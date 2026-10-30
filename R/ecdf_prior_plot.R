# =============================================================================
# ecdf_prior_plot(): compare PRIORS within a single section.
#
# Shares the draws pipeline (.build_ecdf_draws) and renderer (.render_ecdf)
# defined in ecdf_model_plot.R.  Colour = prior; the section is fixed
# (default "All studies").
# =============================================================================

#' ECDF Plot Comparing Priors (fixed section)
#'
#' @description
#' ECDF (+ optional density) comparing the posterior across priors for a single
#' analysis section (default "All studies"). Colour represents the prior.
#' This is the prior-sensitivity companion to \code{\link{ecdf_model_plot}}.
#'
#' @param model A fitted \code{brmsfit} meta-analysis model.
#' @param data Study-level data used to fit \code{model}.
#' @param priors Named list of prior specs; each element a list with
#'   \code{name}, \code{mu_prior}, \code{tau_prior}. The display order of the
#'   legend/colours follows the order of this list (or \code{prior_order}).
#' @param measure Effect measure ("OR","RR","HR","IRR","MD","SMD").
#' @param section_to_plot Which analysis section to display. One of
#'   "All studies", "Excluding High RoB", "PET-PEESE", "Mixture Model", "BMA".
#'   Default "All studies".
#' @param prior_order Optional character vector giving the order of priors in
#'   the legend. Defaults to \code{names(priors)}.
#' @param study_var,rob_var Optional unquoted column names (study id / RoB),
#'   needed only if the chosen \code{section_to_plot} requires them.
#' @param prob_reference "null" or "null_range".
#' @param null_value,null_range,add_null_range,color_null_range Null / ROPE options.
#' @param label_control,label_intervention Arm labels for the density panel.
#' @param title,subtitle,x_label,xlim,x_breaks,color_palette,show_density,font Display options.
#' @param ... Additional section toggles (e.g. \code{exclude_high_rob},
#'   \code{incl_pet_peese}, \code{incl_selection_weight}, \code{incl_robma})
#'   passed on to \code{.build_ecdf_draws()}. Set the one corresponding to
#'   \code{section_to_plot} (other than "All studies") to TRUE so those draws
#'   are built.
#'
#' @return A ggplot (or patchwork) object.
#' @keywords internal
#' @export
#'
#' @seealso \code{\link{ecdf_model_plot}} for comparing model strategies at a fixed prior.
#'
#' @importFrom rlang enquo quo_is_null %||%
#' @importFrom dplyr filter mutate case_when
#' @importFrom purrr map_chr
ecdf_prior_plot <- function(model, data, priors, measure,
                            section_to_plot = "Random Effects",
                            prior_order = NULL,
                            study_var = NULL, rob_var = NULL,
                            prob_reference = "null",
                            null_value = NULL, null_range = NULL,
                            add_null_range = FALSE, color_null_range = "#77bb41",
                            label_control = "Control", label_intervention = "Intervention",
                            title = NULL, subtitle = NULL, x_label = NULL,
                            xlim = NULL, x_breaks = NULL, color_palette = NULL,
                            show_density = TRUE, font = NULL, ...) {

  validate_inputs_sens_plot(model = model, data = data, measure = measure, priors = priors)

  prior_order <- prior_order %||% names(priors)
  if (!all(prior_order %in% names(priors)))
    stop("`prior_order` must contain only names present in `priors`.", call. = FALSE)
  if (!prob_reference %in% c("null", "null_range"))
    stop("`prob_reference` must be one of 'null' or 'null_range'.", call. = FALSE)

  rob_var   <- rlang::enquo(rob_var)
  study_var <- rlang::enquo(study_var)
  props <- get_measure_properties(measure)
  null_value <- null_value %||% props$null_value
  null_range <- .resolve_null_range(measure, null_value, null_range, add_null_range)

  draws <- .build_ecdf_draws(
    model, data, priors, measure,
    rob_var = rob_var, study_var = study_var,
    ...)

  # standardise section labels, then fix to the requested section
  draws <- draws |>
    dplyr::mutate(section_label = dplyr::case_when(
      grepl("^PET-PEESE", section_label) ~ "PET-PEESE",
      grepl("^PET / PEESE", section_label) ~ "PET-PEESE",
      TRUE ~ section_label))

  if (!section_to_plot %in% unique(draws$section_label))
    stop("`section_to_plot` = '", section_to_plot,
         "' not found. Available: ",
         paste(unique(draws$section_label), collapse = ", "),
         ". (Did you enable the matching section toggle?)", call. = FALSE)

  draws_filtered <- dplyr::filter(draws, section_label == section_to_plot)

  # legend/colour order = prior display names in prior_order
  prior_level_order <- unname(purrr::map_chr(priors[prior_order], ~ .x$name %||% .x))

  .render_ecdf(
    draws_filtered, grp_col = "prior_label",
    level_order = prior_level_order, legend_title = "Prior",
    measure = measure, props = props,
    null_value = null_value, null_range = null_range,
    add_null_range = add_null_range, color_null_range = color_null_range,
    prob_reference = prob_reference, color_palette = color_palette,
    xlim = xlim, x_breaks = x_breaks, x_label,
    label_control = label_control, label_intervention = label_intervention,
    title = title, subtitle = subtitle, show_density = show_density, font = font)
}
