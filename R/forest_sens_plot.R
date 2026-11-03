#' Create a Forest Plot with Multi-Prior Pooled-Effect Rows
#'
#' @description
#' A classic forest plot -- individual study rows exactly like
#' [bayes_forest()] -- whose "Pooled Effect" section shows one row per prior
#' specification instead of a single pooled estimate, plus an (optional)
#' prediction-interval row. Useful for showing how sensitive the pooled
#' estimate is to the choice of prior directly on the same plot as the
#' underlying studies. For a fuller sensitivity analysis comparing common-
#' vs random-effects models, RoB exclusion, or RoBMA bias adjustment (but
#' without individual study rows), see [sensitivity_plot()].
#'
#' @param model A fitted Bayesian meta-analysis model object. The study rows
#'   support the same backends as [bayes_forest()]: \code{brmsfit},
#'   \code{meta_stan}, \code{bayesmeta}, \code{stanreg}, \code{stanfit},
#'   \code{CmdStanMCMC}, \code{bayesma}, \code{brma}, \code{rjags}. Refitting
#'   under each entry in \code{priors}, however, only supports
#'   \code{brmsfit}, \code{bayesmeta}, \code{meta_stan} and \code{brma}
#'   (a single \code{RoBMA::brma()} fit, not a \code{RoBMA::RoBMA()}
#'   model-averaged ensemble) -- other backends raise an error naming the
#'   unsupported model class. See \code{\link{bf_priors}} for which prior
#'   families each of those four backends can fit exactly.
#' @param data A data frame containing the study data used for the
#'   meta-analysis.
#' @param measure Character string specifying the effect measure. Must be
#'   one of: "OR", "HR", "RR", "IRR", "MD", or "SMD".
#' @param study_var Column name containing study identifiers/authors.
#' @param year_var Column name containing publication years. Default is NULL.
#' @param n_ctrl Column name containing control group sample sizes.
#' @param n_int Column name containing intervention group sample sizes.
#' @param event_ctrl Column name containing control group event counts.
#' @param event_int Column name containing intervention group event counts.
#' @param mean_ctrl Column name containing control group means.
#' @param mean_int Column name containing intervention group means.
#' @param sd_ctrl Column name containing control group standard deviations.
#' @param sd_int Column name containing intervention group standard
#'   deviations.
#' @param time_ctrl Column name containing control group time periods.
#' @param time_int Column name containing intervention group time periods.
#' @param sort_studies_by Character string specifying how to sort studies.
#'   Options: "author" (default), "year", or "effect".
#' @param priors A named list of prior specifications. Each element should
#'   be a list with \code{name} (display label), \code{mu_prior} (a
#'   \code{bf_prior} object), and optionally \code{tau_prior} (a
#'   \code{bf_prior} object). See \code{\link{bf_normal}} for constructors.
#'   Each prior is refit via the same machinery [sensitivity_plot()] uses,
#'   and shown as its own pooled-effect row labelled "\code{name} Prior".
#' @param prior_order Optional character vector of prior keys specifying
#'   display order. Defaults to the order in \code{priors}.
#' @param primary_prior Optional character key into \code{priors}. That
#'   prior's row label gets " (primary)" appended -- use it to flag which
#'   prior specification matches the prior `model` was actually fit with.
#' @param show_prediction Logical. Show a "Prediction interval" row below
#'   the priors (derived from `model`'s own pooled estimate and tau).
#'   Automatically dropped (with a message) for common-effect models, which
#'   have no tau. Default TRUE.
#' @param study_output Character string specifying how per-study rows are
#'   drawn: "sq_whisker" (default, matching a classic forest plot), "density",
#'   or "boxplot".
#' @param pooled_output Character string specifying how each prior's row is
#'   drawn: "diamond" (default), "density", "boxplot", or "sq_whisker".
#' @param pred_output Character string specifying how the prediction-interval
#'   row is drawn: "diamond" (default), "density", "boxplot", or
#'   "sq_whisker".
#' @param label_control Character string for control group label. Default
#'   "Control".
#' @param label_intervention Character string for intervention group label.
#'   Default "Intervention".
#' @param title Character string for the plot title. Default NULL.
#' @param subtitle Character string for the plot subtitle. Default NULL.
#' @param title_align Character string specifying title alignment.
#' @param xlim Numeric vector of length 2 specifying x-axis limits.
#' @param x_breaks Numeric vector specifying custom x-axis break points.
#' @param x_label Character string for the x-axis label. Default NULL
#'   (measure-specific label).
#' @param null_value Numeric value specifying x-axis value of the null line.
#' @param add_rope Logical indicating whether to add a ROPE. Default FALSE.
#' @param rope_value Numeric vector of length 2 specifying ROPE range.
#' @param rope_color Color for ROPE shading.
#' @param color_palette Character vector of colors for the plot.
#' @param color_study_posterior Color for study posterior densities.
#' @param color_study_posterior_outline Color for study posterior outlines.
#' @param color_pooled_posterior Fill color for each prior's row.
#' @param color_shrinkage_pointinterval Color for shrinkage point intervals.
#' @param color_shrinkage_outline Color for shrinkage plot outlines.
#' @param color_shrinkage_fill Color for shrinkage plot fill.
#' @param split_color_by_null Logical. Split posterior densities by null
#'   value.
#' @param color_favours_control Colour for regions favouring control.
#' @param color_favours_intervention Colour for regions favouring
#'   intervention.
#' @param add_arm_labels Logical. Show "Favours Control/Intervention" labels.
#' @param reverse_arm_labels Logical. Swap the arm label positions.
#' @param color_pred_posterior Color for the prediction-interval row fill.
#' @param color_pred_outline Color for the prediction-interval row outline.
#' @param plot_width Numeric value specifying the relative width of the
#'   plot component. Default 4.
#' @param font Character string specifying the font family.
#'
#' @return A patchwork object containing the complete forest plot.
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
#' forest_sens_plot(
#'   model         = fit,
#'   data          = dat,
#'   measure       = "OR",
#'   study_var      = Author,
#'   n_ctrl        = N_Control,
#'   n_int         = N_Intervention,
#'   event_ctrl    = Event_Control,
#'   event_int     = Event_Intervention,
#'   priors        = priors,
#'   primary_prior = "weakreg"
#' )
#' }
#'
#' @export
forest_sens_plot <- function(model,
                             data,
                             measure,
                             study_var = NULL,
                             year_var = NULL,
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
                             sort_studies_by = "author",
                             priors,
                             prior_order = NULL,
                             primary_prior = NULL,
                             show_prediction = TRUE,
                             study_output = c("density", "boxplot", "sq_whisker"),
                             pooled_output = c("diamond", "density", "boxplot", "sq_whisker"),
                             pred_output = c("diamond", "density", "boxplot", "sq_whisker"),
                             label_control = "Control",
                             label_intervention = "Intervention",
                             title = NULL,
                             subtitle = NULL,
                             title_align = "left",
                             xlim = NULL,
                             x_breaks = NULL,
                             x_label = NULL,
                             null_value = NULL,
                             add_rope = FALSE,
                             rope_value = NULL,
                             rope_color = "grey50",
                             color_palette = NULL,
                             color_study_posterior = "dodgerblue",
                             color_study_posterior_outline = "blue",
                             color_pooled_posterior = "blue",
                             color_shrinkage_pointinterval = "purple",
                             color_shrinkage_outline = "purple",
                             color_shrinkage_fill = NULL,
                             split_color_by_null = FALSE,
                             color_favours_control = "firebrick",
                             color_favours_intervention = "dodgerblue",
                             add_arm_labels = TRUE,
                             reverse_arm_labels = FALSE,
                             color_pred_posterior = "orange",
                             color_pred_outline = "darkorange2",
                             plot_width = 4,
                             font = NULL) {

  validate_inputs(
    model = model, data = data, measure = measure,
    study_var = {{study_var}}, year_var = {{year_var}},
    subgroup = FALSE, subgroup_var = NULL,
    sort_studies_by = sort_studies_by)

  study_output  <- rlang::arg_match(study_output)
  pooled_output <- rlang::arg_match(pooled_output)
  pred_output   <- rlang::arg_match(pred_output)

  # Rename data columns based on measure type (mirrors bayes_forest())
  if (measure %in% c("OR", "RR")) {
    data <- data |>
      dplyr::rename(
        Author = {{study_var}}, Year = {{year_var}},
        N_Control = {{n_ctrl}}, N_Intervention = {{n_int}},
        Event_Control = {{event_ctrl}}, Event_Intervention = {{event_int}})
  } else if (measure %in% c("MD", "SMD")) {
    have_effect <- all(c("yi", "vi") %in% names(data))
    data <- data |>
      dplyr::rename(
        Author = {{study_var}}, Year = {{year_var}},
        N_Control = {{n_ctrl}}, N_Intervention = {{n_int}})
    if (!have_effect) {
      data <- data |>
        dplyr::rename(
          Mean_Control = {{mean_ctrl}}, Mean_Intervention = {{mean_int}},
          SD_Control = {{sd_ctrl}}, SD_Intervention = {{sd_int}})
    }
  } else if (measure == "IRR") {
    data <- data |>
      dplyr::rename(
        Author = {{study_var}}, Year = {{year_var}},
        Time_Control = {{time_ctrl}}, Time_Intervention = {{time_int}},
        Event_Control = {{event_ctrl}}, Event_Intervention = {{event_int}})
  }

  # Disambiguate duplicate Author names (and same-Author/same-Year combos,
  # e.g. multi-arm trials, which otherwise display as identical rows)
  data <- make_authors_unique(data, year = "Year")
  data <- compute_yi_vi(data, measure)

  if (is.null(color_shrinkage_fill)) color_shrinkage_fill <- NA

  model_info  <- detect_model_info(model)
  effect_type <- model_info$pooling
  props       <- get_measure_properties(measure)

  if (!is.null(prior_order)) priors <- priors[prior_order]

  re_flag <- model_has_re(model)
  if (isTRUE(show_prediction) && !isTRUE(re_flag)) {
    cli::cli_inform(c(
      "i" = "{.arg show_prediction} dropped: a common-effect model has no {.field tau},
             so no prediction interval can be derived."))
  }
  show_prediction <- isTRUE(show_prediction) && isTRUE(re_flag)

  pooled_effect <- extract_pooled_summary(model)

  # Study rows only -- the model's own native Pooled Effect row is dropped;
  # it's replaced by the multi-prior block built below.
  study.effect.draws <- forest_data_fn(
    data = data, model = model, subgroup = FALSE,
    sort_studies_by = sort_studies_by, effect_type = effect_type,
    add_pred = FALSE) |>
    dplyr::filter(Author != "Pooled Effect")

  study.summary <- forest_data_summary_fn(
    spread_df = study.effect.draws, data = data, measure = measure,
    sort_studies_by = sort_studies_by, subgroup = FALSE,
    effect_type = effect_type, rob_cols = character())

  study_levels <- levels(droplevels(study.effect.draws$Author))

  # One row per prior, via the same refit machinery sensitivity_plot() uses
  prior_labels <- purrr::imap_chr(priors, function(spec, key) {
    label <- paste0(spec$name %||% key, " Prior")
    if (!is.null(primary_prior) && identical(key, primary_prior)) {
      paste0(label, " (primary)")
    } else {
      label
    }
  })

  prior_draws_raw <- update_prior_sensitivity(
    model = model, data = data, priors = priors, measure = measure,
    section_label = "sensitivity")

  prior_draws <- prior_draws_raw |>
    dplyr::mutate(
      theta_study = if (isTRUE(props$log_scale)) log(x) else x,
      Author = factor(unname(prior_labels[prior]), levels = unname(prior_labels)))

  prior_summary <- prior_draws |>
    dplyr::group_by(Author) |>
    tidybayes::median_qi(theta_study) |>
    dplyr::ungroup() |>
    dplyr::left_join(
      prior_draws |>
        dplyr::group_by(Author) |>
        tidybayes::median_qi(tau) |>
        dplyr::ungroup() |>
        dplyr::select(Author, tau, tau.lower = .lower, tau.upper = .upper),
      by = "Author") |>
    dplyr::transmute(
      Author = as.character(Author),
      Year = NA_integer_,
      weighted_effect = format_sens_estimate(theta_study, .lower, .upper, isTRUE(props$log_scale)),
      unweighted_effect = if (effect_type == "random_effect") {
        paste0("\u03c4 = ", sprintf("%.2f", tau), " [", sprintf("%.2f", tau.lower), ", ", sprintf("%.2f", tau.upper), "]")
      } else {
        weighted_effect
      })

  # Prediction interval, derived from the base model's own pooled estimate + tau
  pred_label <- "Prediction Interval"
  if (isTRUE(show_prediction)) {
    pred_draws <- extract_forest_draws(model, data) |>
      add_prediction_draws() |>
      dplyr::filter(Author == "Prediction") |>
      dplyr::mutate(Author = pred_label)

    pred_summary <- pred_draws |>
      dplyr::group_by(Author) |>
      tidybayes::median_qi(theta_study) |>
      dplyr::ungroup() |>
      dplyr::transmute(
        Author, Year = NA_integer_,
        weighted_effect = format_sens_estimate(theta_study, .lower, .upper, isTRUE(props$log_scale)),
        unweighted_effect = " ")
  } else {
    pred_draws <- NULL
    pred_summary <- NULL
  }

  spacer_label <- " "
  spacer_row <- tibble::tibble(Author = spacer_label, theta_study = NA_real_)
  spacer_summary <- tibble::tibble(
    Author = spacer_label, Year = NA_integer_,
    weighted_effect = " ", unweighted_effect = " ")

  header_row <- tibble::tibble(Author = "Pooled Effect", theta_study = NA_real_)
  header_summary <- tibble::tibble(
    Author = "Pooled Effect", Year = NA_integer_,
    weighted_effect = " ", unweighted_effect = " ")

  author_levels <- c(
    study_levels, spacer_label, "Pooled Effect", unname(prior_labels),
    if (isTRUE(show_prediction)) pred_label)

  combined_draws <- dplyr::bind_rows(
    study.effect.draws, spacer_row, header_row, prior_draws, pred_draws) |>
    dplyr::mutate(Author = factor(as.character(Author), levels = author_levels))

  combined_summary <- dplyr::bind_rows(
    study.summary, spacer_summary, header_summary, prior_summary, pred_summary) |>
    dplyr::mutate(Author = factor(as.character(Author), levels = author_levels))

  plot_break_labels <- c(
    "Pooled Effect",
    if (isTRUE(show_prediction) && length(priors) >= 2) pred_label)

  study.plot <- study.density.plot_fn(
    df = combined_draws,
    pooled_effect = pooled_effect,
    measure = measure,
    effect_type = effect_type,
    subgroup = FALSE,
    color_palette = color_palette,
    color_study_posterior = color_study_posterior,
    color_study_posterior_outline = color_study_posterior_outline,
    color_pooled_posterior = color_pooled_posterior,
    color_overall_posterior = color_pooled_posterior,
    color_shrinkage_outline = color_shrinkage_outline,
    color_shrinkage_pointinterval = color_shrinkage_pointinterval,
    color_shrinkage_fill = color_shrinkage_fill,
    split_color_by_null = split_color_by_null,
    color_favours_control = color_favours_control,
    color_favours_intervention = color_favours_intervention,
    label_control = label_control,
    label_intervention = label_intervention,
    add_arm_labels = add_arm_labels,
    reverse_arm_labels = reverse_arm_labels,
    study_output = study_output,
    pooled_output = pooled_output,
    pooled_labels = unname(prior_labels),
    pred_output = pred_output,
    pred_labels = pred_label,
    add_pred = FALSE,
    color_pred_posterior = color_pred_posterior,
    color_pred_outline = color_pred_outline,
    xlim = xlim,
    x_breaks = x_breaks,
    x_label = x_label,
    null_value = null_value,
    add_rope = add_rope,
    rope_value = rope_value,
    rope_color = rope_color,
    section_break_labels = plot_break_labels,
    font = font)

  table_header_labels <- c("Pooled Effect", if (isTRUE(show_prediction)) pred_label)

  forest.table.left <- forest_table_left_fn(
    forest.data.summary = combined_summary,
    subgroup = FALSE,
    label_control = label_control,
    label_intervention = label_intervention,
    measure = measure,
    pooled_labels = unname(prior_labels),
    header_labels = table_header_labels,
    spacer_labels = "Pooled Effect",
    font = font)

  forest.table.right <- forest_table_right_fn(
    df = combined_summary,
    subgroup = FALSE,
    measure = measure,
    effect_type = effect_type,
    pooled_labels = unname(prior_labels),
    header_labels = table_header_labels,
    spacer_labels = "Pooled Effect",
    font = font)

  patchwork_fn(
    table.left = forest.table.left,
    study.density.plot = study.plot,
    table.right = forest.table.right,
    plot_width = plot_width,
    title = title,
    subtitle = subtitle,
    title_align = title_align,
    font = font)
}

#' Format a point-estimate-with-interval effect string
#'
#' @noRd
format_sens_estimate <- function(theta, lower, upper, log_scale) {
  if (isTRUE(log_scale)) {
    paste0(sprintf("%.2f", exp(theta)), " [", sprintf("%.2f", exp(lower)), ", ", sprintf("%.2f", exp(upper)), "]")
  } else {
    paste0(sprintf("%.2f", theta), " [", sprintf("%.2f", lower), ", ", sprintf("%.2f", upper), "]")
  }
}
