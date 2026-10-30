#' Create posterior plots for Bayesian meta-analysis
#'
#' Produces density plots of the pooled effect (\eqn{\mu}), heterogeneity
#' (\eqn{\tau}), and/or prediction distribution from a fitted model.
#' Supports all bayesfoRest backends: brms, MetaStan, bayesmeta, rstan,
#' cmdstanr, and rstanarm.
#'
#' @param model A fitted model object (brmsfit, meta_stan, bayesmeta,
#'   stanreg, stanfit, or CmdStanMCMC).
#' @param data Data frame with study information. Must include an
#'   \code{Author} column and effect size / SE columns appropriate
#'   for the backend.
#' @param measure Character string for effect measure
#'   ("OR", "RR", "HR", "IRR", "MD", "SMD").
#' @param add_tau Logical, whether to include the heterogeneity (tau) plot.
#' @param add_pred Logical, whether to include the prediction distribution
#'   plot.
#' @param add_summary Logical, whether to annotate the mu plot with a summary
#'   of the effect estimate, 95\% CrI, and posterior probabilities of
#'   benefit/harm. When \code{null_range} is provided, also shows
#'   Pr(Benefit > \eqn{\delta}) and Pr(Harm > \eqn{\delta}). Default
#'   \code{TRUE}.
#' @param color_pred_posterior Color for the prediction posterior
#'   distribution.
#' @param color_pred_posterior_outline Outline color for the prediction
#'   posterior distribution.
#' @param null_value Numeric, the null value for the effect (default from
#'   \code{get_measure_properties}).
#' @param incl_mu_prior Logical, whether to draw the prior for \eqn{\mu} as a
#'   dashed line behind the posterior. Default \code{FALSE}.
#' @param incl_tau_prior Logical, whether to draw the prior for \eqn{\tau} as a
#'   dashed line behind the posterior (needs \code{add_tau = TRUE}). Default
#'   \code{FALSE}.
#' @param mu_prior,tau_prior Optional priors to draw, as \code{bf_prior}
#'   objects (e.g. \code{bf_normal(0, 1)}, \code{bf_half_cauchy(0, 0.5)}) or
#'   density functions on the model scale (log scale for ratio measures). When
#'   \code{NULL}, the priors are read from brms and bayesmeta models. Other
#'   backends do not store their priors, so supply them here.
#' @param color_prior Color for the prior lines.
#' @param null_range Range for the null region: a numeric vector of length 2
#'   (e.g., \code{c(0.9, 1.1)}) or a single half-width around
#'   \code{null_value}. When \code{NULL} and \code{add_null_range = TRUE},
#'   defaults to 0.9 to 1.1 for ratio measures and -0.1 to 0.1 for SMD.
#' @param add_null_range Logical, whether to add shaded null range region.
#' @param color_null_range Color for the null range shading.
#' @param label_control Character string for control group label.
#' @param label_intervention Character string for intervention group label.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param title_align Alignment for title ("left", "center", "right").
#' @param mu_xlim Numeric vector of length 2 for mu plot x-axis limits.
#' @param tau_xlim Numeric vector of length 2 for tau plot x-axis limits.
#' @param pred_xlim Numeric vector of length 2 for prediction plot x-axis
#'   limits.
#' @param x_label Character string for the mu/prediction plot x-axis label.
#' @param x_breaks Numeric vector for x-axis breaks on mu plot.
#' @param tau_breaks Numeric vector for x-axis breaks on tau plot.
#' @param pred_breaks Numeric vector for x-axis breaks on prediction plot.
#' @param color_overall_posterior Color for the overall posterior
#'   distribution.
#' @param color_overall_posterior_outline Outline color for the posterior
#'   distribution.
#' @param split_color_by_null Logical, whether to split posterior color by
#'   null value.
#' @param color_favours_control Color for region favoring control.
#' @param color_favours_intervention Color for region favoring
#'   intervention.
#' @param tau_posterior_color Base color for tau posterior.
#' @param tau_posterior_outline Outline color for tau posterior.
#' @param tau_slab_scale Numeric, scaling factor for tau slab height
#'   relative to mu (default 0.7).
#' @param font Font family for text elements.
#' @param plot_arrangement Character, "vertical" (stacked) or "horizontal"
#'   (side by side).
#' @param ... Additional arguments passed to
#'   \code{extract_forest_draws} (e.g., \code{theta_par},
#'   \code{tau_par}, \code{study_par} for stanfit/CmdStanMCMC models).
#'
#' @details
#' Priors and posteriors are drawn on the same density scale, so a vague prior
#' appears almost flat and a prior that dominates the posterior appears taller
#' than it.
#'
#' @return A ggplot object (or combined plot if \code{add_tau = TRUE} or
#'   \code{add_pred = TRUE}).
#' @export
#'
#' @importFrom grDevices col2rgb rgb
#' @importFrom stats rnorm
#' @importFrom ggplot2 ggplot aes theme_light element_blank element_rect
#'   element_line unit annotate geom_line geom_vline scale_fill_manual
#'   scale_x_log10 scale_x_continuous coord_cartesian labs theme margin
#'   element_text after_stat
#' @importFrom ggdist stat_slab stat_pointinterval scale_fill_ramp_discrete
#'   scale_thickness_shared
#' @importFrom patchwork plot_spacer plot_layout plot_annotation
#' @importFrom dplyr filter group_by summarise mutate n
#' @importFrom tidyr complete
#' @importFrom scales pretty_breaks
#'
overall_plot <- function(model,
                         data,
                         measure,
                         add_tau = FALSE,
                         add_pred = FALSE,
                         add_summary = TRUE,
                         incl_mu_prior = FALSE,
                         incl_tau_prior = FALSE,
                         mu_prior = NULL,
                         tau_prior = NULL,
                         color_prior = "grey40",
                         color_pred_posterior = "orange",
                         color_pred_posterior_outline = "darkorange2",
                         null_value = NULL,
                         null_range = NULL,
                         add_null_range = FALSE,
                         color_null_range = "#77bb41",
                         label_control = "Control",
                         label_intervention = "Intervention",
                         title = NULL,
                         subtitle = NULL,
                         title_align = "left",
                         x_label = NULL,
                         mu_xlim = NULL,
                         tau_xlim = NULL,
                         pred_xlim = NULL,
                         x_breaks = NULL,
                         tau_breaks = NULL,
                         pred_breaks = NULL,
                         color_overall_posterior = "dodgerblue",
                         color_overall_posterior_outline = "blue",
                         split_color_by_null = FALSE,
                         color_favours_control = "firebrick",
                         color_favours_intervention = "dodgerblue",
                         tau_posterior_color = "#77bb41",
                         tau_posterior_outline = NULL,
                         tau_slab_scale = 0.7,
                         font = NULL,
                         plot_arrangement = "vertical",
                         ...) {

  # ========================================================================
  # Validate inputs
  # ========================================================================

  supported <- c("brmsfit", "meta_stan", "bayesmeta", "stanreg",
                 "stanfit", "CmdStanMCMC")
  if (!any(class(model) %in% supported)) {
    cli::cli_abort(
      c("{.arg model} must be one of: {.cls {supported}}.",
        "i" = "Received class: {.cls {class(model)[1]}}")
    )
  }

  if (!"Author" %in% names(data)) {
    cli::cli_abort("{.arg data} must contain an {.field Author} column.")
  }

  # Get measure properties
  props <- get_measure_properties(measure)

  # Set null value from measure properties if not provided
  null_value <- null_value %||% props$null_value
  is_log_scale <- props$log_scale
  x_label <- x_label %||% props$x_label
  null_range <- .resolve_null_range(measure, null_value, null_range,
                                    add_null_range)

  # Set tau x-axis label
  tau_x_label <- if (measure %in% c("MD", "SMD")) {
    "Standard Deviation"
  } else {
    "Standard Deviation (Log Scale)"
  }

  # Set default tau outline color (darker version of base color)
  if (is.null(tau_posterior_outline)) {
    tau_rgb <- grDevices::col2rgb(tau_posterior_color)
    tau_posterior_outline <- grDevices::rgb(
      pmax(0, tau_rgb[1] * 0.5),
      pmax(0, tau_rgb[2] * 0.5),
      pmax(0, tau_rgb[3] * 0.5),
      maxColorValue = 255
    )
  }

  # ========================================================================
  # Extract posterior draws via the backend-agnostic S3 generic
  # ========================================================================

  model_info <- detect_model_info(model, ...)
  is_re <- model_info$pooling == "random_effect"

  # extract_forest_draws returns canonical tibble with:
  #   .draw, Author, theta_study, eta, tau
  draws <- extract_forest_draws(model, data, ...)

  # Filter to pooled effect rows for mu and tau draws

  pooled_draws <- draws |> dplyr::filter(Author == "Pooled Effect")
  mu_samples <- pooled_draws$theta_study

  tau_samples <- if (is_re && !all(is.na(pooled_draws$tau))) {
    pooled_draws$tau
  } else {
    NULL
  }

  # Prior density functions on the model scale (NULL when not drawn)
  model_priors <- if (isTRUE(incl_mu_prior) || isTRUE(incl_tau_prior)) {
    overall_model_priors(model)
  } else {
    list()
  }
  mu_prior_fn <- if (isTRUE(incl_mu_prior)) {
    resolve_prior_density(mu_prior, model_priors$mu, "mu", model)
  }
  tau_prior_fn <- if (isTRUE(incl_tau_prior) && isTRUE(add_tau) &&
                      !is.null(tau_samples)) {
    resolve_prior_density(tau_prior, model_priors$tau, "tau", model)
  }

  # ========================================================================
  # Transform function based on measure type
  # ========================================================================

  transform_fn <- if (is_log_scale) exp else identity

  # --- MU PLOT ---

  # Prepare mu data
  mu_transformed <- transform_fn(mu_samples)
  mu_df <- data.frame(mu = mu_transformed)

  # Set default x limits for mu plot with rounding
  mu_xlim <- if (!is.null(mu_xlim)) {
    mu_xlim
  } else {
    if (is_log_scale) {
      c(0.25, 4)
    } else {
      r <- range(mu_transformed, na.rm = TRUE)
      padding <- 0.1 * diff(r)
      c(floor((r[1] - padding) * 10) / 10,
        ceiling((r[2] + padding) * 10) / 10)
    }
  }

  # Set default x breaks
  if (is.null(x_breaks)) {
    if (is_log_scale) {
      x_breaks <- sort(unique(c(mu_xlim[1], 0.5, 1, 2, mu_xlim[2])))
    } else {
      x_breaks <- scales::pretty_breaks(n = 5)(mu_xlim)
      if (!mu_xlim[1] %in% x_breaks) {
        x_breaks <- sort(unique(c(mu_xlim[1], x_breaks)))
      }
    }
  } else {
    if (!mu_xlim[1] %in% x_breaks) {
      x_breaks <- sort(unique(c(mu_xlim[1], x_breaks)))
    }
  }

  # Determine title alignment
  title_hjust <- switch(title_align,
                        "left" = 0,
                        "center" = 0.5,
                        "right" = 1,
                        0)

  # Define axis line thickness
  axis_line_size <- 0.8

  # Base theme
  base_theme <- ggplot2::theme_light() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white"),
      axis.line.x.bottom = ggplot2::element_line(
        color = "black", linewidth = axis_line_size
      ),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_line(
        color = "black", linewidth = axis_line_size
      ),
      axis.ticks.length = ggplot2::unit(0.15, "cm"),
      plot.title = ggplot2::element_text(hjust = title_hjust),
      plot.subtitle = ggplot2::element_text(hjust = title_hjust),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5, "pt")
    )

  # Start building mu plot
  mu_plot <- ggplot2::ggplot(mu_df, ggplot2::aes(x = mu))

  # Prior line and the slab height that keeps it on the posterior's scale
  mu_overlay <- prior_overlay(mu_prior_fn, mu_samples,
                              lims = if (is_log_scale) log(mu_xlim) else mu_xlim,
                              back = transform_fn, base_scale = 0.9)

  # Add null range shading if requested
  if (add_null_range && !is.null(null_range)) {
    mu_plot <- mu_plot +
      ggplot2::annotate("rect",
                        xmin = null_range[1], xmax = null_range[2],
                        ymin = -Inf, ymax = Inf,
                        fill = color_null_range, alpha = 0.2)
  }

  # Add posterior distribution
  if (split_color_by_null) {
    mu_plot <- mu_plot +
      ggdist::stat_slab(
        ggplot2::aes(fill = ggplot2::after_stat(x < null_value)),
        color = color_overall_posterior_outline,
        alpha = 0.8,
        scale = mu_overlay$slab_scale
      ) +
      ggplot2::scale_fill_manual(
        values = c("TRUE" = color_favours_intervention,
                   "FALSE" = color_favours_control),
        guide = "none"
      )
  } else {
    mu_plot <- mu_plot +
      ggdist::stat_slab(fill = color_overall_posterior,
                        color = color_overall_posterior_outline,
                        alpha = 0.8,
                        scale = mu_overlay$slab_scale)
  }

  # Add point interval
  mu_plot <- mu_plot +
    ggdist::stat_pointinterval(.width = c(.66, .80, .95),
                               color = "black",
                               point_size = 2)

  # Prior on top of the posterior so it stays visible where they overlap
  if (!is.null(mu_overlay$prior_df)) {
    mu_plot <- mu_plot +
      prior_line(mu_overlay$prior_df, color_prior) +
      ggplot2::labs(caption = "Dashed line: prior")
  }

  # Add null reference line (black, solid)
  mu_plot <- mu_plot +
    ggplot2::geom_vline(xintercept = null_value, linetype = "solid",
                        color = "black", linewidth = 0.8)

  # Add direction labels
  label_y <- 0.95
  if (is_log_scale) {
    left_x  <- mu_xlim[1] * 1.4
    right_x <- mu_xlim[2] * 0.65
  } else {
    range_width <- diff(mu_xlim)
    left_x  <- mu_xlim[1] + range_width * 0.12
    right_x <- mu_xlim[2] - range_width * 0.12
  }

  mu_plot <- mu_plot +
    ggplot2::annotate("text", x = left_x, y = label_y,
                      label = paste0("Favours\n", label_intervention),
                      fontface = "bold.italic", size = 3,
                      color = "grey30", vjust = 1) +
    ggplot2::annotate("text", x = right_x, y = label_y,
                      label = paste0("Favours\n", label_control),
                      fontface = "bold.italic", size = 3,
                      color = "grey30", vjust = 1) +
    ggplot2::annotate("segment", x = mu_xlim[1], xend = mu_xlim[2],
                      y = Inf, yend = Inf,
                      linewidth = axis_line_size, color = "grey60")

  # Add summary statistics annotation
  if (isTRUE(add_summary)) {
    # Compute summary on display scale
    est   <- stats::median(mu_transformed, na.rm = TRUE)
    cri   <- stats::quantile(mu_transformed, c(0.025, 0.975), na.rm = TRUE)
    pr_benefit <- round(mean(mu_transformed < null_value) * 100, 1)
    pr_harm    <- round(mean(mu_transformed > null_value) * 100, 1)

    # Format effect estimate and CrI
    fmt <- if (is_log_scale) "%.3f" else "%.2f"
    est_str <- sprintf(fmt, est)
    cri_str <- sprintf(paste0("[", fmt, ", ", fmt, "]"),
                       unname(cri[1]), unname(cri[2]))

    # Build label lines
    summary_lines <- c(
      paste0(measure, " [95% CrI]     ", est_str, "  ", cri_str),
      sprintf("Pr(Benefit)  %5.1f%%       Pr(Harm)  %5.1f%%",
              pr_benefit, pr_harm)
    )

    # Add delta-based probabilities only when null_range is provided
    if (!is.null(null_range)) {
      pr_benefit_delta <- round(
        mean(mu_transformed < null_range[1]) * 100, 1
      )
      pr_harm_delta <- round(
        mean(mu_transformed > null_range[2]) * 100, 1
      )
      summary_lines <- c(
        summary_lines,
        sprintf(
          "Pr(Benefit>\u03b4)  %5.1f%%    Pr(Harm>\u03b4)  %5.1f%%",
          pr_benefit_delta, pr_harm_delta
        )
      )
    }

    summary_label <- paste(summary_lines, collapse = "\n")

    # Position at right edge, mid-height of slab
    summary_x <- if (is_log_scale) mu_xlim[2] * 0.95 else {
      mu_xlim[2] - diff(mu_xlim) * 0.02
    }

    mu_plot <- mu_plot +
      ggplot2::annotate("text",
                        x = summary_x, y = 0.55,
                        label = summary_label,
                        hjust = 1, vjust = 0.5,
                        size = 2.8, color = "grey20",
                        lineheight = 1.2)
  }

  # Scale and coordinate system
  if (is_log_scale) {
    mu_plot <- mu_plot +
      ggplot2::scale_x_log10(breaks = x_breaks,
                             limits = mu_xlim,
                             labels = function(x) sprintf("%.2g", x),
                             expand = c(0, 0)) +
      ggplot2::coord_cartesian(xlim = mu_xlim, ylim = c(0.03, 1),
                               clip = "off")
  } else {
    mu_plot <- mu_plot +
      ggplot2::scale_x_continuous(breaks = x_breaks,
                                  limits = mu_xlim,
                                  labels = function(x) sprintf("%.2g", x),
                                  expand = c(0, 0)) +
      ggplot2::coord_cartesian(xlim = mu_xlim, ylim = c(0.03, 1),
                               clip = "off")
  }

  mu_plot <- mu_plot +
    ggdist::scale_thickness_shared() +
    base_theme +
    ggplot2::labs(x = x_label,
                  y = NULL,
                  title = expression(paste("Overall Effect (", mu, ")")))

  # Apply custom font if specified
  if (!is.null(font)) {
    mu_plot <- mu_plot +
      ggplot2::theme(text = ggplot2::element_text(family = font))
  }

  # ========================================================================
  # PREDICTION PLOT
  # ========================================================================

  pred_plot <- NULL

  if (isTRUE(add_pred)) {
    if (!is_re || is.null(tau_samples)) {
      cli::cli_warn(
        c("Cannot create prediction plot.",
          "i" = "Prediction plots require a random-effects model with tau.")
      )
    } else {
      # Generate prediction draws: theta_pred = mu + N(0, tau)
      pred_draws <- add_prediction_draws(pooled_draws)
      pred_rows  <- pred_draws |> dplyr::filter(Author == "Prediction")
      pred_raw   <- pred_rows$theta_study

      pred_transformed <- transform_fn(pred_raw)
      pred_df <- data.frame(pred = pred_transformed)

      # Set default prediction x limits
      pred_xlim <- pred_xlim %||% mu_xlim

      # Set default prediction breaks
      if (is.null(pred_breaks)) pred_breaks <- x_breaks

      # Build prediction plot
      pred_plot <- ggplot2::ggplot(pred_df, ggplot2::aes(x = pred))

      # Add null range shading if requested
      if (add_null_range && !is.null(null_range)) {
        pred_plot <- pred_plot +
          ggplot2::annotate("rect",
                            xmin = null_range[1], xmax = null_range[2],
                            ymin = -Inf, ymax = Inf,
                            fill = color_null_range, alpha = 0.2)
      }

      # Add posterior distribution
      if (split_color_by_null) {
        pred_plot <- pred_plot +
          ggdist::stat_slab(
            ggplot2::aes(fill = ggplot2::after_stat(x < null_value)),
            color = color_pred_posterior_outline,
            alpha = 0.8
          ) +
          ggplot2::scale_fill_manual(
            values = c("TRUE" = color_favours_intervention,
                       "FALSE" = color_favours_control),
            guide = "none"
          )
      } else {
        pred_plot <- pred_plot +
          ggdist::stat_slab(fill = color_pred_posterior,
                            color = color_pred_posterior_outline,
                            alpha = 0.8)
      }

      # Add point interval
      pred_plot <- pred_plot +
        ggdist::stat_pointinterval(.width = c(.66, .80, .95),
                                   color = "black",
                                   point_size = 2)

      # Add null reference line
      pred_plot <- pred_plot +
        ggplot2::geom_vline(xintercept = null_value, linetype = "solid",
                            color = "black", linewidth = 0.8)

      # Add direction labels
      if (is_log_scale) {
        pred_left_x  <- pred_xlim[1] * 1.4
        pred_right_x <- pred_xlim[2] * 0.65
      } else {
        pred_range_width <- diff(pred_xlim)
        pred_left_x  <- pred_xlim[1] + pred_range_width * 0.12
        pred_right_x <- pred_xlim[2] - pred_range_width * 0.12
      }

      pred_plot <- pred_plot +
        ggplot2::annotate("text", x = pred_left_x, y = label_y,
                          label = paste0("Favours\n", label_intervention),
                          fontface = "bold.italic", size = 3,
                          color = "grey30", vjust = 1) +
        ggplot2::annotate("text", x = pred_right_x, y = label_y,
                          label = paste0("Favours\n", label_control),
                          fontface = "bold.italic", size = 3,
                          color = "grey30", vjust = 1) +
        ggplot2::annotate("segment",
                          x = pred_xlim[1], xend = pred_xlim[2],
                          y = Inf, yend = Inf,
                          linewidth = axis_line_size, color = "grey60")

      # Scale and coordinate system
      if (is_log_scale) {
        pred_plot <- pred_plot +
          ggplot2::scale_x_log10(breaks = pred_breaks,
                                 limits = pred_xlim,
                                 labels = function(x) sprintf("%.2g", x),
                                 expand = c(0, 0)) +
          ggplot2::coord_cartesian(xlim = pred_xlim, ylim = c(0.03, 1),
                                   clip = "off")
      } else {
        pred_plot <- pred_plot +
          ggplot2::scale_x_continuous(breaks = pred_breaks,
                                      limits = pred_xlim,
                                      labels = function(x) sprintf("%.2g", x),
                                      expand = c(0, 0)) +
          ggplot2::coord_cartesian(xlim = pred_xlim, ylim = c(0.03, 1),
                                   clip = "off")
      }

      pred_plot <- pred_plot +
        ggdist::scale_thickness_shared() +
        base_theme +
        ggplot2::labs(x = x_label, y = NULL, title = "Prediction")

      # Apply custom font if specified
      if (!is.null(font)) {
        pred_plot <- pred_plot +
          ggplot2::theme(text = ggplot2::element_text(family = font))
      }
    }
  }

  # ========================================================================
  # TAU PLOT
  # ========================================================================

  tau_plot <- NULL

  if (add_tau) {
    if (is.null(tau_samples)) {
      cli::cli_warn(
        c("No heterogeneity parameter found.",
          "i" = "Tau plot requires a random-effects model.")
      )
    } else {
      tau_df <- data.frame(tau = tau_samples)

      # Define heterogeneity cutpoints
      tau_cuts   <- c(0, 0.1, 0.5, 1, Inf)
      tau_labels_vec <- c("Low", "Reasonable", "Fairly high",
                          "Fairly extreme")

      # Categorize tau samples and calculate percentages
      tau_df$category <- cut(tau_df$tau,
                             breaks = tau_cuts,
                             labels = tau_labels_vec,
                             include.lowest = TRUE,
                             right = FALSE)

      # Calculate percentage in each category
      tau_percentages <- tau_df |>
        dplyr::group_by(category) |>
        dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
        dplyr::mutate(pct = round(n / sum(n) * 100, 1)) |>
        tidyr::complete(
          category = factor(tau_labels_vec, levels = tau_labels_vec),
          fill = list(n = 0, pct = 0)
        )

      pct_low         <- tau_percentages$pct[
        tau_percentages$category == "Low"
      ]
      pct_reasonable  <- tau_percentages$pct[
        tau_percentages$category == "Reasonable"
      ]
      pct_fairly_high <- tau_percentages$pct[
        tau_percentages$category == "Fairly high"
      ]
      pct_extreme     <- tau_percentages$pct[
        tau_percentages$category == "Fairly extreme"
      ]

      # Set default tau x limits with rounding
      tau_xlim <- if (!is.null(tau_xlim)) {
        tau_xlim
      } else {
        r <- range(tau_df$tau, na.rm = TRUE)
        padding <- 0.1 * diff(r)
        c(max(0, floor((r[1] - padding) * 10) / 10),
          ceiling((r[2] + padding) * 10) / 10)
      }

      # Set default tau breaks
      if (is.null(tau_breaks)) {
        reference_points <- c(0, 0.10, 0.25, 0.50, 0.75, 1.00)
        if (tau_xlim[2] > 1.5) {
          extra_breaks <- seq(1.5, tau_xlim[2], by = 0.5)
          extra_breaks <- extra_breaks[extra_breaks < tau_xlim[2]]
          reference_points <- c(reference_points, extra_breaks)
        }
        tau_breaks <- reference_points[reference_points >= tau_xlim[1] &
                                         reference_points <= tau_xlim[2]]
        tau_breaks <- sort(unique(c(tau_xlim[1], tau_breaks, tau_xlim[2])))
      } else {
        if (!tau_xlim[1] %in% tau_breaks)
          tau_breaks <- sort(unique(c(tau_xlim[1], tau_breaks)))
        if (!tau_xlim[2] %in% tau_breaks)
          tau_breaks <- sort(unique(c(tau_breaks, tau_xlim[2])))
      }

      tau_overlay <- prior_overlay(tau_prior_fn, tau_samples,
                                   lims = c(max(tau_xlim[1], 1e-6), tau_xlim[2]),
                                   back = identity, base_scale = tau_slab_scale)

      # Build tau plot with scaled slab height
      tau_plot <- ggplot2::ggplot(tau_df, ggplot2::aes(x = tau))

      tau_plot <- tau_plot +
        ggdist::stat_slab(
          ggplot2::aes(
            fill_ramp = ggplot2::after_stat(
              cut(x,
                  breaks = c(-Inf, 0.1, 0.5, 1, Inf),
                  labels = c("Low", "Reasonable",
                             "Fairly high", "Fairly extreme"))
            )),
          fill = tau_posterior_color,
          color = tau_posterior_outline,
          linewidth = 0.8,
          scale = tau_overlay$slab_scale
        ) +
        ggdist::scale_fill_ramp_discrete(
          range = c(0.8, 0.2),
          guide = "none"
        ) +
        ggdist::stat_pointinterval(.width = c(.66, .80, .95),
                                   color = "black",
                                   point_size = 2)

      if (!is.null(tau_overlay$prior_df)) {
        tau_plot <- tau_plot +
          prior_line(tau_overlay$prior_df, color_prior) +
          ggplot2::labs(caption = "Dashed line: prior")
      }

      # Add heterogeneity category reference lines (dashed, grey)
      tau_plot <- tau_plot +
        ggplot2::geom_vline(xintercept = 0.1, linetype = "dashed",
                            color = "grey40", linewidth = 0.5) +
        ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed",
                            color = "grey40", linewidth = 0.5) +
        ggplot2::geom_vline(xintercept = 1, linetype = "dashed",
                            color = "grey40", linewidth = 0.5)

      very_high_x <- (1 + tau_xlim[2]) / 2

      # Add labels with percentages
      if (tau_xlim[1] < 0.1 && tau_xlim[2] > 0) {
        tau_plot <- tau_plot +
          ggplot2::annotate("text", x = 0.05, y = 0.95,
                            label = paste0("Low\n(", pct_low, "%)"),
                            fontface = "bold.italic", size = 2.6,
                            color = "grey30", vjust = 1)
      }
      if (tau_xlim[1] < 0.5 && tau_xlim[2] > 0.1) {
        tau_plot <- tau_plot +
          ggplot2::annotate("text", x = 0.3, y = 0.95,
                            label = paste0("Moderate\n(",
                                           pct_reasonable, "%)"),
                            fontface = "bold.italic", size = 2.6,
                            vjust = 1)
      }
      if (tau_xlim[1] < 1.0 && tau_xlim[2] > 0.5) {
        tau_plot <- tau_plot +
          ggplot2::annotate("text", x = 0.75, y = 0.95,
                            label = paste0("High\n(",
                                           pct_fairly_high, "%)"),
                            fontface = "bold.italic", size = 2.6,
                            vjust = 1)
      }
      if (tau_xlim[2] > 1.0) {
        tau_plot <- tau_plot +
          ggplot2::annotate("text", x = very_high_x, y = 0.95,
                            label = paste0("Very High\n(",
                                           pct_extreme, "%)"),
                            fontface = "bold.italic", size = 2.6,
                            vjust = 1)
      }

      tau_plot <- tau_plot +
        ggplot2::annotate("segment",
                          x = tau_xlim[1], xend = tau_xlim[2],
                          y = Inf, yend = Inf,
                          linewidth = axis_line_size, color = "grey60")

      # Scale and theme for tau plot
      tau_plot <- tau_plot +
        ggplot2::scale_x_continuous(breaks = tau_breaks,
                                    limits = tau_xlim,
                                    labels = function(x) sprintf("%.2g", x),
                                    expand = c(0, 0)) +
        ggplot2::coord_cartesian(xlim = tau_xlim, ylim = c(0.03, 1),
                                 clip = "off") +
        ggdist::scale_thickness_shared() +
        base_theme +
        ggplot2::labs(
          x = tau_x_label,
          y = NULL,
          title = expression(paste("Heterogeneity (", tau, ")"))
        )

      # Apply custom font if specified
      if (!is.null(font)) {
        tau_plot <- tau_plot +
          ggplot2::theme(text = ggplot2::element_text(family = font))
      }
    }
  }

  # ========================================================================
  # Combine plots
  # ========================================================================

  # Collect all plots that exist
  plot_list <- list(mu_plot)
  if (isTRUE(add_pred) && !is.null(pred_plot)) {
    plot_list <- c(plot_list, list(pred_plot))
  }
  if (isTRUE(add_tau) && !is.null(tau_plot)) {
    plot_list <- c(plot_list, list(tau_plot))
  }

  # If only mu plot, return it directly
  if (length(plot_list) == 1) {
    final_plot <- mu_plot
    if (!is.null(title) || !is.null(subtitle)) {
      final_plot <- final_plot +
        ggplot2::labs(
          title = if (!is.null(title)) title
          else expression(paste("Overall Effect (", mu, ")")),
          subtitle = subtitle
        )
    }
    return(final_plot)
  }

  # Multiple plots -- combine with patchwork
  if (plot_arrangement == "horizontal") {
    spaced_list <- list(plot_list[[1]])
    for (i in seq_along(plot_list)[-1]) {
      spaced_list <- c(spaced_list,
                       list(patchwork::plot_spacer()),
                       list(plot_list[[i]]))
    }
    n_plots   <- length(plot_list)
    n_spacers <- n_plots - 1
    widths <- as.list(rep(c(1, 0.05), n_plots))
    widths <- widths[seq_len(n_plots + n_spacers)]
    combined_plot <- Reduce("+", spaced_list) +
      patchwork::plot_layout(ncol = n_plots + n_spacers,
                             widths = unlist(widths))
  } else {
    combined_plot <- Reduce("/", plot_list) +
      patchwork::plot_layout(ncol = 1)
  }

  # Add title and subtitle using patchwork
  if (!is.null(title) || !is.null(subtitle)) {
    annotation_theme <- ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = title_hjust, face = "bold", size = 14
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = title_hjust, size = 11
      )
    )

    if (!is.null(font)) {
      annotation_theme <- annotation_theme +
        ggplot2::theme(
          text = ggplot2::element_text(family = font),
          plot.title = ggplot2::element_text(
            hjust = title_hjust, face = "bold", size = 14, family = font
          ),
          plot.subtitle = ggplot2::element_text(
            hjust = title_hjust, size = 11, family = font
          )
        )
    }

    combined_plot <- combined_plot +
      patchwork::plot_annotation(
        title    = title,
        subtitle = subtitle,
        theme    = annotation_theme
      )
  }

  combined_plot
}


# =============================================================================
# Internal: prior overlays
# =============================================================================

#' Prior line data and posterior slab height
#'
#' Scales the prior density and the posterior slab by the same factor, so the
#' taller of the two peaks at `base_scale`. Densities are on the model scale;
#' `back` maps model-scale values to the plotted scale.
#'
#' @keywords internal
#' @noRd
prior_overlay <- function(prior_fn, samples, lims, back, base_scale) {
  if (is.null(prior_fn)) {
    return(list(prior_df = NULL, slab_scale = base_scale))
  }
  # ggdist estimates the slab from the draws inside the axis limits
  samples <- samples[samples >= lims[1] & samples <= lims[2]]
  grid <- seq(lims[1], lims[2], length.out = 501)
  prior_y <- prior_fn(grid)
  post_peak <- max(ggdist::density_bounded(samples)$y)
  peak <- max(post_peak, prior_y, na.rm = TRUE)
  list(
    prior_df = data.frame(x = back(grid), y = prior_y / peak * base_scale),
    slab_scale = base_scale * post_peak / peak
  )
}

#' @keywords internal
#' @noRd
prior_line <- function(prior_df, color) {
  ggplot2::geom_line(data = prior_df,
                     ggplot2::aes(x = .data$x, y = .data$y),
                     inherit.aes = FALSE, color = color,
                     linetype = "dashed", linewidth = 0.7)
}

#' Choose the prior to draw: user supplied, else read from the model
#'
#' @keywords internal
#' @noRd
resolve_prior_density <- function(user_prior, model_prior, param, model) {
  if (!is.null(user_prior)) {
    if (is.function(user_prior)) return(user_prior)
    if (inherits(user_prior, "bf_prior")) {
      check_bf_prior(user_prior, param)
      return(bf_prior_density(user_prior))
    }
    cli::cli_abort(
      "{.arg {param}_prior} must be a {.cls bf_prior} (e.g. {.fn bf_normal}) or a density function."
    )
  }
  if (!is.null(model_prior)) return(model_prior)
  cli::cli_warn(c(
    "Could not read the {param} prior from the {.cls {class(model)[1]}} model, so it is not drawn.",
    "i" = "Supply it with {.arg {param}_prior}, e.g. {.code {param}_prior = {if (param == 'mu') 'bf_normal(0, 1)' else 'bf_half_cauchy(0, 0.5)'}}."
  ))
  NULL
}

#' Read the mu and tau priors stored in a fitted model
#'
#' Returns density functions on the model scale. Only brms and bayesmeta store
#' their priors; other backends return NULLs.
#'
#' @keywords internal
#' @noRd
overall_model_priors <- function(model) {
  if (inherits(model, "brmsfit")) {
    pr <- as.data.frame(model$prior)
    pr <- pr[nzchar(pr$prior), , drop = FALSE]
    mu_row <- pr[pr$class == "Intercept" |
                   (pr$class == "b" & pr$coef == "Intercept"), , drop = FALSE]
    tau_row <- pr[pr$class == "sd", , drop = FALSE]
    return(list(
      mu  = if (nrow(mu_row)) brms_prior_density(mu_row$prior[1], half = FALSE),
      tau = if (nrow(tau_row)) brms_prior_density(tau_row$prior[1], half = TRUE)
    ))
  }
  if (inherits(model, "bayesmeta")) {
    mu_pr <- model$mu.prior
    return(list(
      mu  = if (all(is.finite(mu_pr))) {
        function(x) stats::dnorm(x, mu_pr[[1]], mu_pr[[2]])
      },
      tau = if (isTRUE(model$tau.prior.proper)) {
        function(x) model$dprior(tau = x)
      }
    ))
  }
  list(mu = NULL, tau = NULL)
}

#' Density function for a brms prior string such as "student_t(3, 0, 2.5)"
#'
#' @keywords internal
#' @noRd
brms_prior_density <- function(prior_string, half = FALSE) {
  m <- regmatches(prior_string,
                  regexec("^\\s*([a-z_]+)\\((.*)\\)\\s*$", prior_string))[[1]]
  if (length(m) != 3) return(NULL)
  args <- suppressWarnings(as.numeric(strsplit(m[3], ",")[[1]]))
  if (anyNA(args)) return(NULL)
  f <- switch(
    m[2],
    normal      = if (length(args) == 2) function(x) stats::dnorm(x, args[1], args[2]),
    student_t   = if (length(args) == 3) function(x) stats::dt((x - args[2]) / args[3], args[1]) / args[3],
    cauchy      = if (length(args) == 2) function(x) stats::dcauchy(x, args[1], args[2]),
    exponential = if (length(args) == 1) function(x) stats::dexp(x, args[1]),
    lognormal   = if (length(args) == 2) function(x) stats::dlnorm(x, args[1], args[2]),
    uniform     = if (length(args) == 2) function(x) stats::dunif(x, args[1], args[2]),
    NULL
  )
  if (is.null(f) || !half || m[2] %in% c("exponential", "lognormal", "uniform")) {
    return(f)
  }
  truncate_at_zero(f)
}
