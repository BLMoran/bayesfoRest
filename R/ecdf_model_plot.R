# =============================================================================
# ECDF sensitivity plots (brms / bayesfoRest)
#
# Two public functions share one draws pipeline and one renderer:
#   * ecdf_model_plot()  -- fix a PRIOR, compare model strategies/sections
#                           (All studies, Excluding High RoB, PET-PEESE,
#                            Mixture Model, BMA).  Colour = section.
#   * ecdf_prior_plot()  -- fix a SECTION (default "All studies"), compare
#                           PRIORS.  Colour = prior.   [see ecdf_prior_plot.R]
#
# Both build draws via run_ma_sensitivity() (the same backend-agnostic path as
# sensitivity_plot()), so all brms fixes, prior labelling from priors[[]]$name,
# and the exclude_high_rob guard are inherited automatically.
# =============================================================================


# -----------------------------------------------------------------------------
# Internal: build section x prior draws
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.build_ecdf_draws <- function(model, data, priors, measure,
                              rob_var, study_var,
                              exclude_high_rob = FALSE,
                              incl_common_effect = FALSE,
                              incl_random_effect = TRUE,
                              incl_pet_peese = FALSE,
                              incl_selection_weight = FALSE,
                              incl_robma = FALSE,
                              robma_fit = NULL,
                              pet_peese_direction = "negative",
                              pet_peese_threshold = 0.10) {

  prior_labels <- stats::setNames(
    vapply(names(priors), function(nm) priors[[nm]]$name %||% nm, character(1)),
    names(priors))

  sens_specs <- list(
    list(id="common_effect", include = isTRUE(incl_common_effect),
         model=model, data=data, section_label="Common Effect", force_re=FALSE),
    list(id="random_effect", include = isTRUE(incl_random_effect),
         model=model, data=data, section_label="Random Effects", force_re=TRUE),
    list(id="excl_high_rob",
         include = isTRUE(exclude_high_rob) && !rlang::quo_is_null(rob_var),
         model=model,
         data = if (isTRUE(exclude_high_rob) && !rlang::quo_is_null(rob_var))
           dplyr::filter(data, !!rob_var != "High" | is.na(!!rob_var)) else data,
         section_label="Excluding High RoB")
  )

  if (isTRUE(exclude_high_rob) && rlang::quo_is_null(rob_var))
    cli::cli_abort("`exclude_high_rob = TRUE` requires `rob_var`.")
  if (isTRUE(incl_pet_peese) && rlang::quo_is_null(study_var))
    cli::cli_abort("`study_var` must be supplied when `incl_pet_peese = TRUE`.")
  if ((isTRUE(incl_robma) || isTRUE(incl_pet_peese) || isTRUE(incl_selection_weight)) &&
      is.null(robma_fit))
    cli::cli_abort("`robma_fit` must be supplied when any RoBMA option is TRUE.")

  draws_ma <- sens_specs |>
    purrr::keep(~ .x$include) |>
    purrr::map_dfr(run_ma_sensitivity, priors = priors, measure = measure)

  prior_labels <- stats::setNames(
    vapply(names(priors), function(nm) priors[[nm]]$name %||% nm, character(1)),
    names(priors))

  draws_ma <- draws_ma |>
    dplyr::mutate(prior_label = unname(prior_labels[as.character(prior)]))

  # RoBMA-family sections (selection / PET-PEESE / model-averaged), as in sensitivity_plot
  needs_robma <- isTRUE(incl_robma) || isTRUE(incl_pet_peese) || isTRUE(incl_selection_weight)
  draws_robma <- NULL
  if (needs_robma) {
    if (is.null(robma_fit))
      cli::cli_abort("`robma_fit` must be supplied when any RoBMA option is TRUE.")
    robma_types <- c(
      if (isTRUE(incl_pet_peese)) "pet_peese",
      if (isTRUE(incl_selection_weight)) "selection",
      if (isTRUE(incl_robma)) c("conditional","model_averaged")
    )
    robma_fits <- if (inherits(robma_fit, "RoBMA")) list(single = robma_fit) else robma_fit
    draws_robma <- build_robma_sensitivity_draws(
      robma_fits = robma_fits, measure = measure, types = robma_types)
  }

  dplyr::bind_rows(draws_ma, draws_robma)
}


# -----------------------------------------------------------------------------
# Internal: render the ECDF (+ optional density panel)
#   draws_filtered : already reduced to the fixed dimension
#   grp_col        : "section_label" or "prior_label" (the colour aesthetic)
#   level_order    : desired legend/colour order for grp_col
#   legend_title   : legend heading
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.render_ecdf <- function(draws_filtered, grp_col, level_order, legend_title,
                         measure, props,
                         null_value, null_range, add_null_range, color_null_range,
                         prob_reference, color_palette,
                         xlim, x_breaks, x_label,
                         label_control, label_intervention,
                         title, subtitle, show_density, font) {

  # generic grouping column
  draws_filtered$.grp <- draws_filtered[[grp_col]]
  present_levels <- unname(level_order[level_order %in% unique(draws_filtered$.grp)])
  draws_filtered$.grp <- factor(draws_filtered$.grp, levels = present_levels)

  # colours
  if (is.null(color_palette)) {
    n_levels <- length(present_levels)
    n_colors <- max(3, min(9, n_levels))
    color_values <- RColorBrewer::brewer.pal(n_colors, "Set1")[seq_len(n_levels)]
    color_palette <- stats::setNames(unname(color_values), present_levels)
  } else {
    missing_levels <- setdiff(present_levels, names(color_palette))
    if (length(missing_levels) > 0) {
      warning("Some groups missing from color_palette: ",
              paste(missing_levels, collapse = ", "), ". Using default colors.")
      n_missing <- length(missing_levels)
      extra_colors <- RColorBrewer::brewer.pal(max(3, n_missing), "Set2")[seq_len(n_missing)]
      color_palette <- c(color_palette, stats::setNames(extra_colors, missing_levels))
    }
  }

  calc_xlim <- if (!is.null(xlim)) xlim else range(draws_filtered$x, na.rm = TRUE)
  breaks <- x_breaks %||% ggplot2::waiver()
  x_lab <- x_label %||% props$x_label

  if (prob_reference == "null") {
    y_left_label  <- paste0("Probability ", measure, " < ", null_value)
    y_right_label <- paste0("Probability ", measure, " > ", null_value)
  } else {
    y_left_label  <- paste0("Probability ", measure, " < ", null_range[1])
    y_right_label <- paste0("Probability ", measure, " > ", null_range[2])
  }

  ecdf_plot <- ggplot2::ggplot(draws_filtered, ggplot2::aes(x = x, color = .grp)) +
    {if (isTRUE(add_null_range) && !is.null(null_range)) {
      ggplot2::annotate("rect", xmin = null_range[1], xmax = null_range[2],
                        ymin = -Inf, ymax = Inf,
                        fill = scales::alpha(color_null_range, 0.2), color = NA)
    }} +
    {if (isTRUE(add_null_range) && !is.null(null_range)) {
      ggplot2::geom_vline(xintercept = null_range, linetype = "dotted", color = "grey30")
    }} +
    ggplot2::stat_ecdf(geom = "step", linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = null_value, linewidth = 0.8, color = "black") +
    ggplot2::scale_y_continuous(
      name = y_left_label, limits = c(0, 1), breaks = seq(0, 1, 0.1),
      labels = scales::percent_format(),
      sec.axis = ggplot2::sec_axis(~ 1 - ., name = y_right_label,
                                   breaks = seq(0, 1, 0.1),
                                   labels = scales::percent_format())
    ) +
    ggplot2::scale_color_manual(values = color_palette, limits = present_levels,
                                name = legend_title) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = c(0.01, 0.99),
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_rect(fill = "white", color = "grey60", linewidth = 0.3),
      legend.title = ggplot2::element_text(family = font, face = "bold"),
      legend.text = ggplot2::element_text(family = font),
      axis.title = ggplot2::element_text(family = font),
      axis.text = ggplot2::element_text(family = font),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(family = font, face = "bold"),
      plot.subtitle = ggplot2::element_text(family = font),
      axis.line.x.bottom = ggplot2::element_line(linewidth = 0.8, color = "black"),
      axis.ticks.x = ggplot2::element_line(linewidth = 0.8, color = "black"),
      axis.ticks.length.x = ggplot2::unit(0.15, "cm")
    ) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(linetype = "solid", shape = NA)))

  if (isTRUE(props$log_scale)) {
    ecdf_plot <- ecdf_plot +
      ggplot2::scale_x_log10(breaks = breaks, expand = c(0, 0)) +
      ggplot2::coord_cartesian(xlim = calc_xlim, clip = "off") +
      ggplot2::labs(x = x_lab)
  } else {
    ecdf_plot <- ecdf_plot +
      ggplot2::scale_x_continuous(breaks = breaks, expand = c(0, 0)) +
      ggplot2::coord_cartesian(xlim = calc_xlim, clip = "off") +
      ggplot2::labs(x = x_lab)
  }

  if (!isTRUE(show_density)) return(ecdf_plot)

  ecdf_plot <- ecdf_plot +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x.bottom = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = 2, l = 5.5, unit = "pt")
    )

  density_plot <- ggplot2::ggplot(draws_filtered, ggplot2::aes(x = x, color = .grp)) +
    {if (isTRUE(add_null_range) && !is.null(null_range)) {
      ggplot2::annotate("rect", xmin = null_range[1], xmax = null_range[2],
                        ymin = -Inf, ymax = Inf,
                        fill = scales::alpha(color_null_range, 0.2), color = NA)
    }} +
    {if (isTRUE(add_null_range) && !is.null(null_range)) {
      ggplot2::geom_vline(xintercept = null_range, linetype = "dotted", color = "grey30")
    }} +
    ggdist::stat_slab(fill = NA, slab_linewidth = 0.9, normalize = "all") +
    ggplot2::geom_vline(xintercept = null_value, linewidth = 0.8, color = "black") +
    ggplot2::scale_color_manual(values = color_palette, limits = present_levels,
                                name = legend_title) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(family = font),
      axis.text.x = ggplot2::element_text(family = font),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line.x.bottom = ggplot2::element_line(linewidth = 0.8, color = "black"),
      axis.ticks.x = ggplot2::element_line(linewidth = 0.8, color = "black"),
      axis.ticks.length.x = ggplot2::unit(0.15, "cm"),
      plot.margin = ggplot2::margin(t = 2, r = 5.5, b = 5.5, l = 5.5, unit = "pt")
    ) +
    ggplot2::labs(x = x_lab)

  if (isTRUE(props$log_scale)) {
    density_plot <- density_plot +
      ggplot2::scale_x_log10(breaks = breaks, expand = c(0, 0)) +
      ggplot2::coord_cartesian(xlim = calc_xlim, clip = "off")
  } else {
    density_plot <- density_plot +
      ggplot2::scale_x_continuous(breaks = breaks, expand = c(0, 0)) +
      ggplot2::coord_cartesian(xlim = calc_xlim, clip = "off")
  }

  density_plot <- density_plot +
    ggplot2::annotation_custom(
      grid::textGrob(
        label = paste0("Favours\n", label_intervention),
        x = grid::unit(0.05, "npc"), y = grid::unit(0.6, "npc"),
        just = c("left", "bottom"),
        gp = grid::gpar(col = "grey30", fontsize = 9, fontface = "bold.italic", fontfamily = font)
      ), xmin = calc_xlim[1] - 0.01, xmax = calc_xlim[2], ymin = -Inf, ymax = Inf
    ) +
    ggplot2::annotation_custom(
      grid::textGrob(
        label = paste0("Favours\n", label_control),
        x = grid::unit(0.97, "npc"), y = grid::unit(0.6, "npc"),
        just = c("right", "bottom"),
        gp = grid::gpar(col = "grey30", fontsize = 9, fontface = "bold.italic", fontfamily = font)
      ), xmin = calc_xlim[1], xmax = calc_xlim[2] - 0.01, ymin = -Inf, ymax = Inf
    )

  (ecdf_plot / patchwork::plot_spacer() / density_plot) +
    patchwork::plot_layout(ncol = 1, heights = c(1, 0.0, 0.35),
                           axes = "collect_x", axis_titles = "collect_x") &
    ggplot2::coord_cartesian(xlim = calc_xlim)
}


# -----------------------------------------------------------------------------
# Internal: shared null-range resolution
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.resolve_null_range <- function(measure, null_value, null_range, add_null_range) {
  if (is.null(null_range) && isTRUE(add_null_range)) {
    null_range <- switch(
      measure,
      OR = c(0.9, 1.1), RR = c(0.9, 1.1), HR = c(0.9, 1.1), IRR = c(0.9, 1.1),
      SMD = c(-0.1, 0.1),
      stop("For MD, `null_range` must be supplied.", call. = FALSE)
    )
  } else if (!is.null(null_range)) {
    null_range <- if (length(null_range) == 1)
      c(null_value - null_range, null_value + null_range) else sort(null_range)
  }
  null_range
}


# =============================================================================
# PUBLIC: compare model strategies / sections at a fixed prior
# =============================================================================
#' ECDF Plot Comparing Model Strategies (fixed prior)
#'
#' @description
#' ECDF (+ optional density) comparing sensitivity model strategies -- All
#' studies, Excluding High RoB, PET-PEESE, Mixture Model, BMA -- for a single
#' prior. Colour represents the model strategy/section.
#'
#' @param model A fitted \code{brmsfit} meta-analysis model.
#' @param data Study-level data used to fit \code{model}.
#' @param priors Named list of prior specs; each element a list with
#'   \code{name}, \code{mu_prior}, \code{tau_prior}.
#' @param measure Effect measure ("OR","RR","HR","IRR","MD","SMD").
#' @param prior_to_plot Name (in \code{priors}) of the prior to display.
#'   Defaults to the first prior.
#' @param study_var,rob_var Optional unquoted column names (study id / RoB).
#' @param exclude_high_rob,incl_common_effect,incl_random_effect,incl_pet_peese,incl_selection_weight,incl_robma
#'   Logical section toggles.
#' @param robma_fit A fitted RoBMA object (required if \code{incl_robma}).
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
#' @importFrom dplyr filter mutate bind_rows case_when select any_of left_join
#' @importFrom purrr keep map_dfr imap_dfr
#' @importFrom ggplot2 ggplot aes annotate geom_vline stat_ecdf scale_y_continuous sec_axis scale_color_manual theme_light theme element_rect element_text element_blank element_line unit margin labs guides guide_legend scale_x_log10 scale_x_continuous coord_cartesian expansion waiver annotation_custom
#' @importFrom ggdist stat_slab
#' @importFrom scales alpha percent_format
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats setNames
#' @importFrom grid textGrob unit gpar
#' @importFrom patchwork plot_spacer plot_layout
ecdf_model_plot <- function(model,
                            data,
                            priors,
                            measure,
                            prior_to_plot = NULL,
                            study_var = NULL,
                            rob_var = NULL,
                            exclude_high_rob = FALSE,
                            incl_common_effect = FALSE,
                            incl_random_effect = TRUE,
                            incl_pet_peese = FALSE,
                            incl_selection_weight = FALSE,
                            incl_robma = FALSE,
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

  validate_inputs_sens_plot(model = model, data = data, measure = measure, priors = priors)

  prior_to_plot <- prior_to_plot %||% names(priors)[[1]]
  if (!prior_to_plot %in% names(priors))
    stop("`prior_to_plot` must be one of: ", paste(names(priors), collapse = ", "), call. = FALSE)
  if (!prob_reference %in% c("null", "null_range"))
    stop("`prob_reference` must be one of 'null' or 'null_range'.", call. = FALSE)

  rob_var   <- rlang::enquo(rob_var)
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
    robma_fit,
    pet_peese_direction,
    pet_peese_threshold)

  draws_filtered <- draws |>
    dplyr::filter(prior == prior_to_plot) |>
    dplyr::mutate(section_label = dplyr::case_when(
      grepl("^PET-PEESE", section_label) ~ "PET-PEESE",
      grepl("^PET / PEESE", section_label) ~ "PET-PEESE",
      TRUE ~ section_label))

  section_order <- c("Common Effect", "Random Effects", "Excluding High RoB",
                     "PET-PEESE", "Selection Model", "RoBMA (conditional)",
                     "RoBMA (model-averaged)")

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
