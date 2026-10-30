#' Generate Sensitivity Analysis Plot for Bayesian Meta-Analysis
#'
#' @description
#' Creates a comprehensive sensitivity analysis visualization showing how
#' meta-analytic estimates vary across different prior specifications,
#' analysis approaches, and publication bias adjustments. Supports all
#' model backends (brms, MetaStan, bayesmeta, rstan, cmdstanr) and
#' integrates RoBMA for bias-adjusted sensitivity analyses.
#'
#' @param model A fitted model object (brmsfit, meta_stan, bayesmeta,
#'   stanreg, stanfit, or CmdStanMCMC).
#' @param data A data frame containing the study data.
#' @param priors A named list of prior specifications. Each element should
#'   be a list with \code{name} (display label), \code{mu_prior} (a
#'   \code{bf_prior} object), and optionally \code{tau_prior} (a
#'   \code{bf_prior} object). See \code{\link{bf_normal}} for constructors.
#' @param measure Character. Effect measure: "OR", "RR", "HR", "IRR",
#'   "MD", or "SMD".
#' @param prior_order Optional character vector of prior keys specifying display
#'   order. Defaults to the order in \code{priors}.
#' @param model_order Optional character vector of section labels specifying
#'   display order. Overrides the default ordering.
#' @param rob_var Optional. Name of risk of bias variable (unquoted).
#' @param exclude_high_rob Logical. Include sensitivity excluding high RoB
#'   studies. Default FALSE.
#' @param incl_common_effect Logical. Include common-effect model section.
#'   Default FALSE. Only applies when the base model is random-effects.
#' @param incl_random_effect Logical. Include random-effects model section.
#'   Default TRUE.
#' @param incl_selection_weight Logical. Include RoBMA selection (weight)
#'   model rows. Default FALSE. Requires \code{robma_fit}.
#' @param incl_pet_peese Logical. Include RoBMA PET-PEESE rows. Default FALSE.
#'   Requires \code{robma_fit}.
#' @param incl_robma Logical. Include full RoBMA sensitivity rows (conditional
#'   and model-averaged). Default FALSE. Requires \code{robma_fit}.
#' @param robma_fit A single fitted RoBMA object (from \code{RoBMA::RoBMA()})
#'   or a named list of RoBMA objects (one per prior
#'   level). Required when any RoBMA option is TRUE.
#' @param seed Integer seed for reproducibility. Default TRUE (uses model seed).
#' @param add_probs Logical. Add probability columns. Default FALSE.
#' @param null_value Numeric. Null value for hypothesis testing.
#' @param null_range Numeric vector of length 2 or single value. ROPE.
#' @param add_null_range Logical. Show null range region. Default FALSE.
#' @param color_null_range Character. Color for null range. Default "#77bb41".
#' @param label_control Character. Control group label. Default "Control".
#' @param label_intervention Character. Intervention label. Default
#'   "Intervention".
#' @param title Optional title.
#' @param subtitle Optional subtitle.
#' @param title_align Character. "left", "center", or "right".
#' @param xlim Optional x-axis limits.
#' @param x_breaks Optional x-axis breaks.
#' @param color_palette Optional named vector of colors for sections.
#' @param color_overall_posterior Fill color. Default "dodgerblue".
#' @param color_overall_posterior_outline Outline color. Default "blue".
#' @param split_color_by_null Logical. Color by null. Default FALSE.
#' @param color_favours_control Character. Default "firebrick".
#' @param color_favours_intervention Character. Default "dodgerblue".
#' @param plot_width Numeric. Width ratio for density plot. Default 4.
#' @param font Optional font family.
#'
#' @return A patchwork object combining tables and density plots.
#'
#' @examples
#' \dontrun{
#' priors <- list(
#'   vague = list(
#'     name      = "Vague",
#'     mu_prior  = bf_normal(0, 10),
#'     tau_prior = bf_half_cauchy(0, 1)
#'   ),
#'   weakreg = list(
#'     name      = "Weakly Regularising",
#'     mu_prior  = bf_normal(0, 1),
#'     tau_prior = bf_half_cauchy(0, 0.5)
#'   ),
#'   informative = list(
#'     name      = "Informative",
#'     mu_prior  = bf_normal(0, 0.5),
#'     tau_prior = bf_half_cauchy(0, 0.3)
#'   )
#' )
#'
#' sensitivity_plot(
#'   model   = fit,
#'   data    = dat,
#'   priors  = priors,
#'   measure = "OR",
#'   incl_common_effect = TRUE,
#'   exclude_high_rob   = TRUE,
#'   rob_var            = rob_rating
#' )
#' }
#'
#' @keywords internal
#' @export
sensitivity_plot <- function(model,
                             data,
                             priors,
                             measure,
                             prior_order = NULL,
                             model_order = NULL,
                             rob_var = NULL,
                             exclude_high_rob = FALSE,
                             incl_common_effect = FALSE,
                             incl_random_effect = TRUE,
                             incl_selection_weight = FALSE,
                             incl_pet_peese = FALSE,
                             incl_robma = FALSE,
                             robma_fit = NULL,
                             seed = TRUE,
                             add_probs = FALSE,
                             null_value = NULL,
                             null_range = NULL,
                             add_null_range = FALSE,
                             color_null_range = "#77bb41",
                             label_control = "Control",
                             label_intervention = "Intervention",
                             title = NULL,
                             subtitle = NULL,
                             title_align = "left",
                             xlim = NULL,
                             x_breaks = NULL,
                             color_palette = NULL,
                             color_overall_posterior = "dodgerblue",
                             color_overall_posterior_outline = "blue",
                             split_color_by_null = FALSE,
                             color_favours_control = "firebrick",
                             color_favours_intervention = "dodgerblue",
                             plot_width = 4,
                             font = NULL) {

  # ===========================
  # 1. Validation
  # ===========================
  validate_inputs_sens_plot(
    model   = model,
    data    = data,
    measure = measure,
    priors  = priors
  )

  # Check RoBMA availability
  needs_robma <- isTRUE(incl_robma) || isTRUE(incl_pet_peese) || isTRUE(incl_selection_weight)
  if (needs_robma && is.null(robma_fit)) {
    stop(
      "`robma_fit` must be supplied when any RoBMA option is TRUE ",
      "(incl_robma, incl_pet_peese, incl_selection_weight).",
      call. = FALSE
    )
  }

  rob_var <- rlang::enquo(rob_var)
  if (isTRUE(exclude_high_rob) && rlang::quo_is_null(rob_var)) {
    cli::cli_abort("`exclude_high_rob = TRUE` requires `rob_var` to name the RoB column.")
  }

  props      <- get_measure_properties(measure)
  null_value <- null_value %||% props$null_value

  # ===========================
  # 2. Null range handling
  # ===========================
  if (is.null(null_range) && isTRUE(add_null_range)) {
    null_range <- switch(
      measure,
      OR = c(0.9, 1.1), RR = c(0.9, 1.1), HR = c(0.9, 1.1), IRR = c(0.9, 1.1),
      SMD = c(-0.1, 0.1),
      stop("For MD, `null_range` must be supplied.", call. = FALSE)
    )
  } else if (!is.null(null_range)) {
    null_range <- if (length(null_range) == 1) {
      c(null_value - null_range, null_value + null_range)
    } else {
      sort(null_range)
    }
  }

  # ===========================
  # 3. Prior order
  # ===========================
  if (!is.null(prior_order)) {
    priors <- priors[prior_order]
  }

  # ===========================
  # 4. Build MA draws (standard sections)
  # ===========================
  model_info <- detect_model_info(model)

  sens_specs <- list(
    # Common-effect section
    list(
      id            = "common_effect",
      include       = isTRUE(incl_common_effect),
      model         = model,
      data          = data,
      section_label = "Common Effect",
      force_re      = FALSE
    ),

    # Random-effects section
    list(
      id            = "random_effect",
      include       = isTRUE(incl_random_effect),
      model         = model,
      data          = data,
      section_label = "Random Effects",
      force_re      = TRUE
    ),

    # Excluding high RoB
    list(
      id            = "excl_high_rob",
      include = isTRUE(exclude_high_rob),
      data    = if (isTRUE(exclude_high_rob))
        dplyr::filter(data, !!rob_var != "High" | is.na(!!rob_var)) else data,
      model         = model,
      section_label = "Excluding High RoB"
    ))

  draws_ma <- sens_specs |>
    purrr::keep(~ .x$include) |>
    purrr::map_dfr(run_ma_sensitivity,
                   priors  = priors,
                   measure = measure)

  # ===========================
  # 5. Build RoBMA draws
  # ===========================
  draws_robma <- NULL

  if (needs_robma && !is.null(robma_fit)) {

    # Normalise robma_fit
    is_single_robma <- inherits(robma_fit, "RoBMA")
    if (is_single_robma) {
      # Single fit: wrap with a single key — draws appear once per section
      robma_fits <- list(single = robma_fit)
    } else {
      robma_fits <- robma_fit
    }

    # Determine which RoBMA types to include
    robma_types <- character()
    if (isTRUE(incl_pet_peese))        robma_types <- c(robma_types, "pet_peese")
    if (isTRUE(incl_selection_weight)) robma_types <- c(robma_types, "selection")
    if (isTRUE(incl_robma)) {
      robma_types <- c(robma_types, "conditional", "model_averaged")
    }

    section_labels <- c(
      pet_peese      = "PET-PEESE",
      selection      = "Selection Model (Weight)",
      conditional    = "RoBMA (Conditional)",
      model_averaged = "RoBMA (Model-Averaged)"
    )

    draws_robma <- build_robma_sensitivity_draws(
      robma_fits     = robma_fits,
      measure        = measure,
      types          = robma_types,
      section_labels = section_labels
    )
  }

  draws <- dplyr::bind_rows(draws_ma, draws_robma)

  # ===========================
  # 6. Apply model_order
  # ===========================
  if (!is.null(model_order)) {
    draws <- draws |>
      dplyr::mutate(
        section_label = factor(section_label, levels = model_order)
      ) |>
      dplyr::filter(!is.na(section_label))
  }

  # ===========================
  # 7. Prior table
  # ===========================
  prior_table_ma <- extract_mu_tau_priors(priors = priors, model = model)

  prior_table_robma <- tibble::tibble(
    prior = character(), prior_label = character(),
    mu_prior_unicode = character(), tau_prior_unicode = character()
  )

  if (!is.null(draws_robma) && nrow(draws_robma) > 0) {
    prior_table_robma <- purrr::imap_dfr(robma_fits, function(fit, prior_name) {
      tryCatch({
        priors_info <- extract_priors_from_robma_fit(fit)
        tibble::tibble(
          prior             = prior_name,
          prior_label       = if (is_single_robma) "" else (priors[[prior_name]]$name %||% prior_name),
          mu_prior_unicode  = priors_info$mu_prior_unicode,
          tau_prior_unicode = priors_info$tau_prior_unicode
        )
      }, error = function(e) {
        tibble::tibble(
          prior = prior_name,
          prior_label = if (is_single_robma) "" else (priors[[prior_name]]$name %||% prior_name),
          mu_prior_unicode = NA_character_,
          tau_prior_unicode = NA_character_
        )
      })
    })
  }

  prior_table <- dplyr::bind_rows(prior_table_ma, prior_table_robma)

  # ===========================
  # 8. Posterior summary
  # ===========================
  summary_df <- draws |>
    summarise_sensitivity_posteriors(
      null_value = null_value,
      null_range = null_range
    )

  # Determine which sections are RoBMA vs MA for prior matching
  robma_sections <- if (!is.null(draws_robma) && nrow(draws_robma) > 0) {
    unique(draws_robma$section_label)
  } else {
    character(0)
  }

  summary_df <- summary_df |>
    dplyr::mutate(
      mu_prior_unicode = dplyr::if_else(
        section_label %in% robma_sections,
        prior_table_robma$mu_prior_unicode[
          match(paste(prior, prior_label),
                paste(prior_table_robma$prior, prior_table_robma$prior_label))
        ],
        prior_table_ma$mu_prior_unicode[
          match(paste(prior, prior_label),
                paste(prior_table_ma$prior, prior_table_ma$prior_label))
        ]
      ),
      tau_prior_unicode = dplyr::if_else(
        section_label %in% robma_sections,
        prior_table_robma$tau_prior_unicode[
          match(paste(prior, prior_label),
                paste(prior_table_robma$prior, prior_table_robma$prior_label))
        ],
        prior_table_ma$tau_prior_unicode[
          match(paste(prior, prior_label),
                paste(prior_table_ma$prior, prior_table_ma$prior_label))
        ]
      )
    ) |>
    dplyr::mutate(
      tau_prior_unicode = dplyr::if_else(!is.na(re) & re, tau_prior_unicode, NA_character_)
    )

  # ===========================
  # 9. Tables
  # ===========================
  table_left <- sensitivity_table_left(summary_df, font = font)
  table_right <- sensitivity_table_right(
    summary_df,
    measure    = measure,
    add_probs  = add_probs,
    null_range = null_range,
    font       = font
  )

  # ===========================
  # 10. Density plot
  # ===========================
  sens_plot <- sensitivity_density_plot_fn(
    df                             = draws,
    measure                        = measure,
    split_color_by_null            = split_color_by_null,
    color_overall_posterior         = color_overall_posterior,
    color_overall_posterior_outline = color_overall_posterior_outline,
    color_favours_control          = color_favours_control,
    color_favours_intervention     = color_favours_intervention,
    label_control                  = label_control,
    label_intervention             = label_intervention,
    xlim                           = xlim,
    x_breaks                       = x_breaks,
    null_value                     = null_value,
    null_range                     = null_range,
    color_null_range               = color_null_range,
    add_null_range                 = add_null_range,
    font                           = font
  )

  # ===========================
  # 11. Patchwork assembly
  # ===========================
  patchwork_fn(
    table.left         = table_left,
    study.density.plot = sens_plot,
    table.right        = table_right,
    plot_width         = plot_width,
    title              = title,
    subtitle           = subtitle,
    title_align        = title_align,
    font               = font
  )
}
