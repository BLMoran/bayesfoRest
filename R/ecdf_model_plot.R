#' ECDF Plot Comparing MA Models ( with a fixed prior)
#'
#' @description
#' ECDF (+ optional density) comparing sensitivity model strategies -- All
#' studies, Excluding High RoB, PET-PEESE, Mixture Model, BMA -- for a single
#' prior. Colour represents the model strategy/section.
#'
#' @param model A fitted \code{brmsfit} meta-analysis model.
#' @param data Study-level data used to fit \code{model}.
#' @param prior A single prior spec: a list with \code{name}, \code{mu_prior},
#'   and optionally \code{tau_prior}. To compare across multiple priors
#'   instead of model strategies, use \code{\link{ecdf_prior_plot}}.
#' @param measure Effect measure ("OR","RR","HR","IRR","MD","SMD").
#' @param study_var,rob_var Optional unquoted column names (study id / RoB).
#' @param exclude_high_rob,incl_common_effect,incl_random_effect,incl_pet_peese,incl_selection_weight,incl_robma,incl_robma_bma
#'   Logical section toggles. \code{incl_robma} adds the RoBMA conditional
#'   row (pooled effect given an effect exists); \code{incl_robma_bma} adds
#'   the model-averaged row (can show a point mass at the null) -- the two
#'   are independent.
#' @param robma_fit A single fitted RoBMA object for this prior (required if
#'   \code{incl_robma}, \code{incl_robma_bma}, \code{incl_pet_peese}, or
#'   \code{incl_selection_weight}).
#' @param pet_peese_direction,pet_peese_threshold PET-PEESE options.
#' @param prob_reference "null" or "null_range".
#' @param null_value,null_range,add_null_range,color_null_range Null / ROPE options.
#' @param label_control,label_intervention Arm labels for the density panel.
#' @param title,subtitle,x_label,xlim,x_breaks,color_palette,show_density,font Display options.
#'
#' @return A ggplot (or patchwork) object.
#' @keywords internal
#' @export
#'
#' @importFrom rlang enquo quo_is_null %||%
#' @importFrom dplyr filter mutate bind_rows case_when select any_of left_join group_by group_modify ungroup across all_of
#' @importFrom purrr keep map_dfr imap_dfr
#' @importFrom ggplot2 ggplot aes annotate geom_vline geom_line scale_y_continuous sec_axis scale_color_manual theme_light theme element_rect element_text element_blank element_line unit margin labs guides guide_legend scale_x_log10 scale_x_continuous coord_cartesian expansion waiver annotation_custom
#' @importFrom ggdist stat_slab
#' @importFrom scales alpha percent_format
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats setNames
#' @importFrom grid textGrob unit gpar
#' @importFrom patchwork plot_spacer plot_layout
ecdf_model_plot <- function(model,
                            data,
                            prior,
                            measure,
                            study_var = NULL,
                            rob_var = NULL,
                            exclude_high_rob = FALSE,
                            incl_common_effect = FALSE,
                            incl_random_effect = TRUE,
                            incl_pet_peese = FALSE,
                            incl_selection_weight = FALSE,
                            incl_robma = FALSE,
                            incl_robma_bma = FALSE,
                            robma_fit = NULL,
                            pet_peese_direction = "negative",
                            pet_peese_threshold = 0.10,
                            prob_reference = "null",
                            null_value = NULL,
                            null_range = NULL,
                            add_null_range = FALSE,
                            color_null_range = "#77bb41",
                            label_control = "Control",
                            label_intervention = "Intervention",
                            title = NULL,
                            subtitle = NULL,
                            x_label = NULL,
                            xlim = NULL,
                            x_breaks = NULL,
                            color_palette = NULL,
                            show_density = TRUE,
                            font = NULL) {

  priors <- stats::setNames(list(prior), "single")
  validate_inputs_sens_plot(model = model, data = data, measure = measure, priors = priors)

  if (!prob_reference %in% c("null", "null_range"))
    stop("`prob_reference` must be one of 'null' or 'null_range'.", call. = FALSE)

  needs_robma <- isTRUE(incl_robma) || isTRUE(incl_robma_bma) ||
    isTRUE(incl_pet_peese) || isTRUE(incl_selection_weight)
  if (needs_robma && !is.null(robma_fit) && !inherits(robma_fit, "RoBMA")) {
    cli::cli_abort(c(
      "{.arg robma_fit} must be a single fitted RoBMA object for {.fn ecdf_model_plot}.",
      "i" = "Pass the fit for this prior directly, e.g. {.code robma_fit = robma_fit_weakreg}.",
      "i" = "To compare RoBMA results across multiple priors, use {.fn ecdf_prior_plot} instead."))
  }

  rob_var <- rlang::enquo(rob_var)
  study_var <- rlang::enquo(study_var)
  props <- get_measure_properties(measure)
  null_value <- null_value %||% props$null_value
  null_range <- .resolve_null_range(measure, null_value, null_range, add_null_range)

  draws <- .build_ecdf_draws(
    model, data, priors, measure,
    rob_var, study_var,
    exclude_high_rob,
    incl_common_effect,
    incl_random_effect,
    incl_pet_peese,
    incl_selection_weight,
    incl_robma,
    incl_robma_bma,
    robma_fit,
    pet_peese_direction,
    pet_peese_threshold)

  draws_filtered <- draws |>
    dplyr::mutate(section_label = dplyr::case_when(
      grepl("^PET-PEESE", section_label) ~ "PET-PEESE",
      grepl("^PET / PEESE", section_label) ~ "PET-PEESE",
      TRUE ~ section_label))

  section_order <- c("Common Effect", "Random Effects", "Excluding High RoB",
                     "PET-PEESE", "Selection Model (Weight)", "RoBMA (Conditional)",
                     "RoBMA (Model-Averaged)")

  .render_ecdf(
    draws_filtered, grp_col = "section_label",
    level_order = section_order, legend_title = "Analysis",
    measure = measure, props = props,
    null_value = null_value, null_range = null_range,
    add_null_range = add_null_range, color_null_range = color_null_range,
    prob_reference = prob_reference, color_palette = color_palette,
    xlim = xlim, x_breaks = x_breaks, x_label,
    label_control = label_control, label_intervention = label_intervention,
    title = title, subtitle = subtitle, show_density = show_density, font = font)
}
