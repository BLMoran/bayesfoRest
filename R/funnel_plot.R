#' Funnel Plot (multi-backend)
#'
#' @description
#' Publication-ready funnel plot for meta-analysis, displaying study-level
#' effect sizes against their precision (standard error), styled to match the
#' bayesfoRest package. The plot itself only needs study-level \code{yi}/\code{sei}
#' and a single pooled estimate, so it accepts any fitted model backend that
#' bayesfoRest supports for the pooled line: \code{brmsfit}, \code{bayesmeta},
#' MetaStan / \code{stanfit} (rstan), \code{CmdStanMCMC} (cmdstanr),
#' \code{rjags}/\code{runjags}, or \code{bayesma}. The pooled value is obtained
#' via the same backend-agnostic accessor used by \code{bayes_forest()}.
#'
#' @param model A fitted meta-analysis model (any supported backend). Used only
#'   to draw the pooled-estimate line (when \code{pooled_line = "posterior"}).
#' @param data A data frame of study-level data. Must either contain the columns
#'   named by \code{yi_col}/\code{sei_col}, or the arm columns needed to compute
#'   \code{yi}/\code{sei} for \code{measure} (as in \code{bayes_forest()}).
#' @param measure Character; effect measure. One of \code{"OR"}, \code{"RR"},
#'   \code{"HR"}, \code{"IRR"}, \code{"MD"}, \code{"SMD"}.
#' @param studyvar Unquoted column with study identifiers.
#' @param yearvar Unquoted column with study year (optional, used in labels).
#' @param yi_col,sei_col Optional character names of pre-computed effect and SE
#'   columns (e.g. difference-of-medians \code{yi}/\code{sei}). If supplied,
#'   \code{yi}/\code{sei} are taken directly and no arm columns are needed.
#' @param n_ctrl,n_int,event_ctrl,event_int,mean_ctrl,mean_int,sd_ctrl,sd_int,time_ctrl,time_int
#'   Unquoted arm columns, used only when \code{yi_col}/\code{sei_col} are not
#'   supplied.
#' @param subgroup Logical; colour points by subgroup.
#' @param subgroupvar Unquoted subgroup column (when \code{subgroup = TRUE}).
#' @param contour_lines Logical; add 90/95/99% significance contour shading.
#' @param contour_alpha Numeric; base alpha for contour shading.
#' @param pooled_line One of \code{"posterior"} (posterior median of the pooled
#'   effect, backend-agnostic), \code{"fixed"} (inverse-variance fixed-effect
#'   estimate from the data), or \code{"none"}.
#' @param color_points,color_points_outline,color_pooled,color_null,color_contour
#'   Colours.
#' @param point_size,label_size Sizes.
#' @param xlim,x_breaks,x_label Axis limits / breaks / label.
#' @param null_value,add_null_line Null-effect line.
#' @param label_studies Logical; label points with study (year).
#' @param title,subtitle,title_align,font Display options.
#'
#' @return A ggplot object.
#' @export
#'
#' @importFrom rlang enquo quo_is_null eval_tidy arg_match %||%
#' @importFrom dplyr mutate filter if_else rename bind_rows
#' @importFrom purrr map list_rbind
#' @importFrom stats qnorm median weighted.mean
funnel_plot <- function(model,
                        data,
                        measure,
                        studyvar = NULL,
                        yearvar = NULL,
                        yi_col = NULL,
                        sei_col = NULL,
                        n_ctrl = NULL,
                        n_int = NULL,
                        event_ctrl = NULL,
                        event_int = NULL,
                        mean_ctrl = NULL,
                        mean_int = NULL,
                        sd_ctrl = NULL,
                        sd_int = NULL,
                        time_ctrl = NULL,
                        time_int = NULL,
                        subgroup = FALSE,
                        subgroupvar = NULL,
                        contour_lines = TRUE,
                        contour_alpha = 0.15,
                        pooled_line = c("posterior", "fixed", "none"),
                        color_points = "dodgerblue",
                        color_points_outline = "blue",
                        color_pooled = "blue",
                        color_null = "black",
                        color_contour = "grey50",
                        point_size = 3,
                        xlim = NULL,
                        x_breaks = NULL,
                        null_value = NULL,
                        add_null_line = FALSE,
                        label_studies = FALSE,
                        label_size = 3,
                        title = NULL,
                        subtitle = NULL,
                        x_label = NULL,
                        title_align = "left",
                        font = NULL) {

  pooled_line <- rlang::arg_match(pooled_line)

  props      <- get_measure_properties(measure)
  x_lab <- x_label %||% props$x_label
  null_value <- null_value %||% props$null_value

  studyvar_q     <- rlang::enquo(studyvar)
  year_q         <- rlang::enquo(yearvar)
  subgroup_var_q <- rlang::enquo(subgroupvar)

  # ---------------------------------------------------------------------------
  # 1. Study-level yi / sei
  #    Either taken directly (yi_col/sei_col) or computed from arm columns.
  # ---------------------------------------------------------------------------
  if (!is.null(yi_col) && !is.null(sei_col)) {
    data$yi  <- data[[yi_col]]
    data$sei <- data[[sei_col]]
  } else {
    es <- data |>
      dplyr::rename(
        N_Control = {{n_ctrl}},
        N_Intervention = {{n_int}},
        Event_Control = {{event_ctrl}},
        Event_Intervention = {{event_int}},
        Mean_Control = {{mean_ctrl}},
        Mean_Intervention = {{mean_int}},
        SD_Control = {{sd_ctrl}},
        SD_Intervention = {{sd_int}},
        Time_Control = {{time_ctrl}},
        Time_Intervention = {{time_int}}
      ) |>
      compute_yi_vi(measure)
    data$yi  <- es$yi
    data$sei <- sqrt(es$vi)
  }

  # study labels / year / subgroup, resolved without assuming column names
  data$Author <- if (!rlang::quo_is_null(studyvar_q))
    as.character(rlang::eval_tidy(studyvar_q, data)) else as.character(data$Author)
  data$Year <- if (!rlang::quo_is_null(year_q))
    rlang::eval_tidy(year_q, data) else if ("Year" %in% names(data)) data$Year else NA
  if (isTRUE(subgroup) && !rlang::quo_is_null(subgroup_var_q))
    data$Subgroup <- as.character(rlang::eval_tidy(subgroup_var_q, data))

  # ---------------------------------------------------------------------------
  # 2. Pooled estimate (backend-agnostic)
  # ---------------------------------------------------------------------------
  pooled_estimate <- switch(
    pooled_line,
    none      = NA_real_,
    fixed     = stats::weighted.mean(data$yi, 1 / data$sei^2, na.rm = TRUE),
    # funnel_pooled() returns the outcome scale; the plot works on the model
    # (log) scale for ratio measures, like the "fixed" estimate above
    posterior = {
      pooled <- funnel_pooled(model, data, measure)
      if (isTRUE(props$log_scale)) log(pooled) else pooled
    }
  )

  # ---------------------------------------------------------------------------
  # 3. Plot (unchanged from the single-backend version)
  # ---------------------------------------------------------------------------
  plot_data <- data |>
    dplyr::mutate(
      x_val = if (isTRUE(props$log_scale)) exp(yi) else yi,
      Author_label = dplyr::if_else(
        is.na(Year),
        as.character(Author),
        paste0(Author, " (", Year, ")")
      )
    )

  # ---- Axis limits ----
  if (isTRUE(props$log_scale)) {
    if (is.null(xlim)) {
      x_min <- min(exp(data$yi - 2 * data$sei), na.rm = TRUE) * 0.8
      x_max <- max(exp(data$yi + 2 * data$sei), na.rm = TRUE) * 1.2
      xlim <- c(max(0.05, x_min), x_max)
    }
  } else {
    if (is.null(xlim)) {
      x_range <- max(abs(data$yi) + 2 * data$sei, na.rm = TRUE)
      xlim <- c(-x_range, x_range)
    }
  }

  breaks <- if (!is.null(x_breaks)) x_breaks else ggplot2::waiver()
  se_max <- max(data$sei, na.rm = TRUE) * 1.1

  # ---- Build funnel contours (centred on pooled estimate) ----
  contour_data <- NULL
  z_vals    <- c(qnorm(0.995), qnorm(0.975), qnorm(0.95))  # 2.576, 1.96, 1.645
  ci_labels <- c("99% CI", "95% CI", "90% CI")
  se_seq    <- seq(0, se_max, length.out = 300)

  if (isTRUE(contour_lines) && !is.na(pooled_estimate)) {
    contour_data <- purrr::map(seq_along(z_vals), function(i) {
      if (isTRUE(props$log_scale)) {
        lower <- exp(pooled_estimate - z_vals[i] * se_seq)
        upper <- exp(pooled_estimate + z_vals[i] * se_seq)
      } else {
        lower <- pooled_estimate - z_vals[i] * se_seq
        upper <- pooled_estimate + z_vals[i] * se_seq
      }
      tibble::tibble(
        x = c(lower, rev(upper)),
        y = c(se_seq, rev(se_seq)),
        level = ci_labels[i]
      )
    }) |>
      purrr::list_rbind() |>
      dplyr::mutate(level = factor(level, levels = ci_labels))
  }

  # ---- Build the plot ----
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_val, y = sei))

  if (isTRUE(contour_lines) && !is.null(contour_data)) {
    alpha_steps <- c(contour_alpha, contour_alpha * 1.5, contour_alpha * 2.5)

    for (i in seq_along(ci_labels)) {
      cd <- contour_data |> dplyr::filter(level == ci_labels[i])
      p <- p +
        ggplot2::geom_polygon(
          data = cd,
          ggplot2::aes(x = x, y = y),
          fill = color_contour,
          alpha = alpha_steps[i],
          color = NA,
          inherit.aes = FALSE
        )
    }

    border_data <- purrr::map(seq_along(z_vals), function(i) {
      if (isTRUE(props$log_scale)) {
        lower <- exp(pooled_estimate - z_vals[i] * se_seq)
        upper <- exp(pooled_estimate + z_vals[i] * se_seq)
      } else {
        lower <- pooled_estimate - z_vals[i] * se_seq
        upper <- pooled_estimate + z_vals[i] * se_seq
      }
      dplyr::bind_rows(
        tibble::tibble(x = lower, y = se_seq, side = "lower", level = ci_labels[i]),
        tibble::tibble(x = upper, y = se_seq, side = "upper", level = ci_labels[i])
      )
    }) |> purrr::list_rbind()

    p <- p +
      ggplot2::geom_line(
        data = border_data,
        ggplot2::aes(x = x, y = y, group = interaction(level, side)),
        color = color_contour, linetype = "dashed", linewidth = 0.4,
        inherit.aes = FALSE
      )
  }

  # Null line
  if (isTRUE(add_null_line)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = null_value,
        color = color_null, linewidth = 0.75
      )
  }

  # Pooled estimate line
  if (pooled_line != "none" && !is.na(pooled_estimate)) {
    pooled_x <- if (isTRUE(props$log_scale)) exp(pooled_estimate) else pooled_estimate
    p <- p +
      ggplot2::geom_vline(
        xintercept = pooled_x,
        color = color_pooled, linewidth = 0.75, linetype = "solid"
      )
  }

  # Study points
  if (isTRUE(subgroup) && "Subgroup" %in% names(plot_data)) {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(fill = Subgroup),
        shape = 21, size = point_size,
        color = color_points_outline, stroke = 0.6
      ) +
      ggplot2::scale_fill_brewer(palette = "Set2", name = "Subgroup")
  } else {
    p <- p +
      ggplot2::geom_point(
        fill = color_points, color = color_points_outline,
        shape = 21, size = point_size, stroke = 0.6
      )
  }

  # Study labels
  if (isTRUE(label_studies)) {
    rlang::check_installed("ggrepel", reason = "to label studies.")
    p <- p +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = Author_label),
        size = label_size, color = "grey30",
        max.overlaps = 20, family = font
      )
  }

  # ---- Reverse y-axis (SE: 0 at top, max at bottom) ----
  p <- p +
    ggplot2::scale_y_reverse(
      expand = c(0, 0),
      limits = c(se_max, 0)
    )

  # ---- X-axis scale ----
  if (isTRUE(props$log_scale)) {
    p <- p +
      ggplot2::scale_x_log10(
        breaks = breaks,
        expand = c(0, 0),
        limits = xlim
      )
  } else {
    p <- p +
      ggplot2::scale_x_continuous(
        breaks = breaks,
        expand = c(0, 0),
        limits = xlim
      )
  }

  # ---- Theme ----
  p <- p +
    ggplot2::theme_light() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey92"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line.x.bottom = ggplot2::element_line(color = "black", linewidth = 0.75),
      axis.line.y.left = ggplot2::element_line(color = "black", linewidth = 0.75),
      axis.text = ggplot2::element_text(color = "black", family = font),
      axis.title = ggplot2::element_text(family = font),
      plot.margin = ggplot2::margin(10, 10, 10, 10),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold", family = font),
      legend.text = ggplot2::element_text(family = font)
    ) +
    ggplot2::labs(
      x = x_lab,,
      y = "Standard Error"
    )

  # ---- Title / subtitle ----
  if (!is.null(title) || !is.null(subtitle)) {
    hjust_val <- switch(title_align,
                        "left" = 0,
                        "center" = 0.5,
                        "centre" = 0.5,
                        "right" = 1,
                        0)

    if (!is.null(title)) {
      p <- p + ggplot2::labs(title = title) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = 16, face = "bold", hjust = hjust_val,
            margin = ggplot2::margin(b = if (is.null(subtitle)) 10 else 5),
            family = font
          )
        )
    }
    if (!is.null(subtitle)) {
      p <- p + ggplot2::labs(subtitle = subtitle) +
        ggplot2::theme(
          plot.subtitle = ggplot2::element_text(
            size = 14, hjust = hjust_val,
            margin = ggplot2::margin(b = 10),
            color = "gray30", family = font
          )
        )
    }
  }

  p
}


# =============================================================================
# Backend-agnostic pooled estimate for the funnel line.
# Mirrors bayes_forest()'s model handling: one generic, per-backend methods,
# and a default that reuses extract_pooled_draws() (so any backend already
# wired for bayes_forest works here for free).
# =============================================================================

#' @keywords internal
#' @noRd
funnel_pooled <- function(model, data, measure) UseMethod("funnel_pooled")

# ratio measures live on the log scale internally -> exponentiate the median
.funnel_ratio <- function(measure) measure %in% c("OR", "RR", "HR", "IRR")

#' @keywords internal
#' @noRd
.funnel_median <- function(mu, measure) {
  stats::median(if (.funnel_ratio(measure)) exp(mu) else mu)
}

#' @export
#' @noRd
funnel_pooled.brmsfit <- function(model, data, measure) {
  mu <- posterior::as_draws_df(model)$b_Intercept
  .funnel_median(mu, measure)
}

#' @export
#' @noRd
funnel_pooled.bayesmeta <- function(model, data, measure) {
  mu <- model$summary["median", "mu"]
  if (.funnel_ratio(measure)) exp(mu) else mu
}

# rstan / MetaStan (stanfit) and cmdstanr (CmdStanMCMC) share the draws pattern
#' @keywords internal
#' @noRd
.funnel_pooled_stan <- function(draws_df, measure,
                                par = c("mu", "pooled", "b_Intercept", "theta_mu")) {
  nm <- intersect(par, names(draws_df))
  if (length(nm) == 0)
    cli::cli_abort("Could not find a pooled-effect parameter in the model draws.
                    Looked for: {paste(par, collapse = ', ')}.")
  .funnel_median(draws_df[[nm[1]]], measure)
}

#' @export
#' @noRd
funnel_pooled.stanfit <- function(model, data, measure)
  .funnel_pooled_stan(posterior::as_draws_df(model), measure)

#' @export
#' @noRd
funnel_pooled.CmdStanMCMC <- function(model, data, measure)
  .funnel_pooled_stan(posterior::as_draws_df(model$draws()), measure)

#' @export
#' @noRd
funnel_pooled.stanreg <- function(model, data, measure)
  .funnel_pooled_stan(posterior::as_draws_df(model), measure)

# rjags / runjags (coda mcmc.list)
#' @keywords internal
#' @noRd
.funnel_pooled_jags <- function(mcmc, measure,
                                par = c("mu", "pooled", "theta.mu", "d")) {
  m  <- as.matrix(mcmc)
  nm <- intersect(par, colnames(m))
  if (length(nm) == 0)
    cli::cli_abort("Could not find a pooled-effect node in the JAGS output.
                    Looked for: {paste(par, collapse = ', ')}.")
  .funnel_median(m[, nm[1]], measure)
}

#' @export
#' @noRd
funnel_pooled.runjags <- function(model, data, measure) {
  rlang::check_installed("coda")
  .funnel_pooled_jags(coda::as.mcmc.list(model), measure)
}

#' @export
#' @noRd
funnel_pooled.mcmc.list <- function(model, data, measure)
  .funnel_pooled_jags(model, measure)

#' @export
#' @noRd
funnel_pooled.jags <- function(model, data, measure) {
  rlang::check_installed("coda")
  .funnel_pooled_jags(coda::as.mcmc.list(model$mcmc %||% model), measure)
}

#' @export
#' @noRd
funnel_pooled.bayesma <- function(model, data, measure) {
  mu <- extract_pooled_summary(model)[["estimate"]]
  if (.funnel_ratio(measure)) exp(mu) else mu
}

# Fallback: reuse the same extractor bayes_forest()/sensitivity_plot() dispatch on.
# Any backend already supported there works here with no extra method.
#' @export
#' @noRd
funnel_pooled.default <- function(model, data, measure) {
  stats::median(extract_pooled_draws(model, data, measure))
}
