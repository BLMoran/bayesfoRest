#' Create ECDF Plot for Sensitivity Analysis
#'
#' @description
#' Creates an empirical cumulative distribution function (ECDF) plot showing 
#' posterior distributions across different sensitivity analyses (All studies, 
#' Excluding High RoB, PET-PEESE, Mixture Model, BMA).
#'
#' @param model A fitted brmsfit object from a Bayesian meta-analysis
#' @param data A data frame containing the study data used to fit the model
#' @param priors A named list of brms prior specifications (created with \code{brms::prior()})
#'   with names typically including "vague", "weakreg", and "informative"
#' @param measure Character string specifying the effect measure.  One of "OR" (odds ratio),
#'   "RR" (risk ratio), "HR" (hazard ratio), "IRR" (incidence rate ratio), 
#'   "MD" (mean difference), or "SMD" (standardized mean difference)
#' @param study_var Optional.  Name of the study identifier variable (unquoted). Required
#'   if \code{incl_pet_peese = TRUE} or \code{incl_mixture = TRUE}
#' @param rob_var Optional. Name of the risk of bias variable (unquoted). Should contain
#'   values like "Low", "Some concerns", "High"
#' @param exclude_high_rob Logical. If TRUE, performs sensitivity analysis excluding
#'   studies with high risk of bias.  Default is FALSE
#' @param incl_pet_peese Logical. If TRUE, includes PET-PEESE bias adjustment analysis.
#'   Default is FALSE
#' @param pet_peese_direction Character string. Direction of expected bias for PET-PEESE.  
#'   Either "negative" or "positive".  Default is "negative"
#' @param pet_peese_threshold Numeric. Threshold for choosing between PET and PEESE 
#'   (based on proportion of studies showing effect). Default is 0.10
#' @param incl_mixture Logical. If TRUE, includes robust mixture model analysis.
#'   Default is FALSE
#' @param incl_bma Logical. If TRUE, includes Bayesian model averaging results.
#'   Default is FALSE
#' @param model_bma Named list of RoBMA model fits (from \code{RoBMA_brms()}) to include
#'   in model averaging. Names should be "vague", "weakreg", or "informative"
#' @param prior_to_plot Character string. Which prior specification to plot in the ECDF. 
#'   One of "vague", "weakreg", or "informative". Default is "weakreg"
#' @param prob_reference Character string. What reference to use for probability calculations.
#'   Either "null" (compare to null_value) or "null_range" (compare to null_range boundaries).
#'   Default is "null"
#' @param null_value Numeric. The null value for the effect measure. If NULL, uses
#'   default based on measure (1 for ratio measures, 0 for difference measures)
#' @param null_range Numeric vector of length 2 or single numeric value.  Defines the
#'   range of practical equivalence around the null value. 
#' @param add_null_range Logical. If TRUE, displays the null range region on the plot.
#'   Default is FALSE
#' @param color_null_range Character string. Color for the null range region.
#'   Default is "#77bb41" (green)
#' @param label_control Character string. Label for the control/reference group.
#'   Default is "Control"
#' @param label_intervention Character string. Label for the intervention group.
#'   Default is "Intervention"
#' @param title Optional character string. Main title for the plot
#' @param subtitle Optional character string. Subtitle for the plot
#' @param xlim Optional numeric vector of length 2.  X-axis limits for the plot
#' @param x_breaks Optional numeric vector.  Custom x-axis break points
#' @param color_palette Optional named character vector of colors for each section.  
#'   Names should match section labels (e.g., "All studies", "Excluding High RoB", etc.)
#' @param show_density Logical. If TRUE, includes a density plot below the ECDF. 
#'   Default is TRUE
#' @param font Optional character string. Font family to use for the plot
#'
#' @return A ggplot object (or patchwork object if \code{show_density = TRUE})
#'
#' @details
#' This function creates an ECDF plot that displays the cumulative probability 
#' distributions of pooled effects across different sensitivity analyses.  The plot
#' includes:
#' \itemize{
#'   \item Left y-axis: Probability that mean effect < x
#'   \item Right y-axis: Probability that mean effect > x
#'   \item Optional density curves below the ECDF
#'   \item Optional region of practical equivalence (ROPE) shading
#' }
#' 
#' The \code{prob_reference} argument controls how probabilities are calculated:
#' \itemize{
#'   \item "null":  Probabilities are calculated relative to the null value (e.g., OR = 1)
#'   \item "null_range":  Probabilities are calculated relative to the null range boundaries
#'     (left y-axis shows P(effect < lower bound), right y-axis shows P(effect > upper bound))
#' }
#'
#' @examples
#' \dontrun{
#' # Basic ECDF plot
#' sensitivity_ecdf_plot(
#'   model = my_brms_model,
#'   data = my_data,
#'   priors = list(
#'     vague = prior("normal(0, 10)", class = "Intercept"),
#'     weakreg = prior("normal(0, 1)", class = "Intercept"),
#'     informative = prior("normal(0.5, 0.5)", class = "Intercept")
#'   ),
#'   measure = "OR",
#'   prior_to_plot = "weakreg"
#' )
#'
#' # Advanced ECDF with null range probability reference
#' sensitivity_ecdf_plot(
#'   model = my_model,
#'   data = my_data,
#'   priors = my_priors,
#'   measure = "OR",
#'   study_var = study_id,
#'   rob_var = rob,
#'   exclude_high_rob = TRUE,
#'   incl_pet_peese = TRUE,
#'   incl_mixture = TRUE,
#'   incl_bma = TRUE,
#'   model_bma = my_robma_models,
#'   prior_to_plot = "weakreg",
#'   prob_reference = "null_range",
#'   null_range = c(0.9, 1.1),
#'   add_null_range = TRUE
#' )
#' }
#'
#' @export
#' 
ecdf_plot <- function(model,
                      data,
                      priors,
                      measure,
                      study_var = NULL,
                      rob_var = NULL,
                      exclude_high_rob = FALSE,
                      incl_pet_peese = FALSE,
                      pet_peese_direction = "negative",
                      pet_peese_threshold = 0.10,
                      incl_mixture = FALSE,
                      incl_bma = FALSE,
                      model_bma = NULL,
                      prob_reference = "null",
                      null_value = NULL,
                      null_range = NULL,
                      add_null_range = FALSE,
                      color_null_range = "#77bb41",
                      label_control = "Control",
                      label_intervention = "Intervention",
                      title = NULL,
                      subtitle = NULL,
                      xlim = NULL,
                      x_breaks = NULL,
                      color_palette = NULL,
                      show_density = TRUE,
                      prior_to_plot = "weakreg",
                      font = NULL) {
  
  # ---------------------------
  # 1. Validation
  # ---------------------------
  validate_inputs_sens_plot(
    model   = model,
    data    = data,
    measure = measure,
    priors  = priors
  )
  
  if (incl_bma && is.null(model_bma)) {
    stop("`model_bma` must be supplied when `incl_bma = TRUE`.", call. = FALSE)
  }
  
  if (! prob_reference %in% c("null", "null_range")) {
    stop("`prob_reference` must be one of 'null' or 'null_range'.", call. = FALSE)
  }
  
  # Capture tidy eval variables
  rob_var <- rlang::enquo(rob_var)
  study_var <- rlang::enquo(study_var)
  
  # Validate study_var if PET-PEESE or mixture is included
  if (isTRUE(incl_pet_peese) && rlang::quo_is_null(study_var)) {
    stop("`study_var` must be supplied when `incl_pet_peese = TRUE`.", call. = FALSE)
  }
  
  if (isTRUE(incl_mixture) && rlang::quo_is_null(study_var)) {
    stop("`study_var` must be supplied when `incl_mixture = TRUE`.", call. = FALSE)
  }
  
  props <- get_measure_properties(measure)
  null_value <- null_value %||% props$null_value
  
  # ---------------------------
  # 2. Null range handling
  # ---------------------------
  if (is.null(null_range) && isTRUE(add_null_range)) {
    null_range <- switch(
      measure,
      OR  = c(0.9, 1.1),
      RR  = c(0.9, 1.1),
      HR  = c(0.9, 1.1),
      IRR = c(0.9, 1.1),
      SMD = c(-0.1, 0.1),
      stop("For MD, `null_range` must be supplied.", call. = FALSE)
    )
  } else if (! is.null(null_range)) {
    null_range <- if (length(null_range) == 1) {
      c(null_value - null_range, null_value + null_range)
    } else {
      sort(null_range)
    }
  }
  
  # ---------------------------
  # 3. Prior labels
  # ---------------------------
  prior_labels <- c(
    vague       = "Vague",
    weakreg     = "Weakly Regularising",
    informative = "Informative"
  )
  
  # ---------------------------
  # 4. Define section order for legend
  # ---------------------------
  section_order <- c(
    "All studies",
    "Excluding High RoB",
    "PET-PEESE",
    "Mixture Model",
    "BMA"
  )
  
  # ---------------------------
  # 5. Build draws (ALL SECTIONS)
  # ---------------------------
  sens_specs <- list(
    list(
      id = "all",
      include = TRUE,
      model = model,
      data = data,
      section_label = "All studies"
    ),
    list(
      id = "excl_high_rob",
      include = isTRUE(exclude_high_rob),
      model = model,
      data = dplyr::filter(data, !! rob_var != "High" | is.na(!!rob_var)),
      section_label = "Excluding High RoB"
    ),
    list(
      id = "pet_peese",
      include = isTRUE(incl_pet_peese),
      model = model,
      data = data,
      section_label = "PET-PEESE",
      study_var = study_var,
      pet_peese_direction = pet_peese_direction,
      pet_peese_threshold = pet_peese_threshold
    ),
    list(
      id = "mixture",
      include = isTRUE(incl_mixture),
      model = model,
      data = data,
      section_label = "Mixture Model",
      study_var = study_var
    )
  )
  
  draws_ma <- sens_specs |>
    purrr::keep(~ .x$include) |>
    purrr::map_dfr(run_ma_sensitivity,
                   priors = priors,
                   measure = measure,
                   prior_labels = prior_labels)
  
  if (isTRUE(incl_bma)) {
    draws_bma <- purrr::imap_dfr(
      model_bma,
      ~ robma_to_sensitivity_draws(
        robma_fit   = .x,
        prior       = .y,
        prior_label = c(
          vague = "Vague",
          weakreg = "Weakly Regularising",
          informative = "Informative"
        )[[.y]],
        section_label = "BMA"))
  } else {
    draws_bma <- NULL
  }
  
  draws <- dplyr::bind_rows(draws_ma, draws_bma)

  # ---------------------------
  # 6. Filter to selected prior
  # ---------------------------
  draws_filtered <- draws |>
    dplyr::filter(prior == prior_to_plot)
  
  # ---------------------------
  # 7. Get actual sections present in the data and standardize labels
  # ---------------------------
  actual_sections <- unique(draws_filtered$section_label)
  
  # Standardize section labels - handle PET-PEESE variations
  draws_filtered <- draws_filtered |>
    dplyr::mutate(
      section_label = dplyr::case_when(
        grepl("^PET-PEESE", section_label) ~ "PET-PEESE",
        grepl("^PET / PEESE", section_label) ~ "PET-PEESE",
        TRUE ~ section_label
      )
    )
  
  # Update actual_sections after standardization
  actual_sections <- unique(draws_filtered$section_label)
  
  # ---------------------------
  # 8. Set section_label as ordered factor
  # ---------------------------
  # Only include sections that are present in the data, maintaining order
  present_sections <- section_order[section_order %in% actual_sections]
  
  draws_filtered <- draws_filtered |>
    dplyr::mutate(
      section_label = factor(section_label, levels = present_sections)
    )

  # ---------------------------
  # 9. Set up colors - ensure all present sections have colors
  # ---------------------------
  if (is.null(color_palette)) {
    # Use RColorBrewer Set1 palette like the example
    n_sections <- length(present_sections)
    # Set1 has 9 colors max, but we need at least 3 for brewer.pal
    n_colors <- max(3, min(9, n_sections))
    color_values <- RColorBrewer::brewer.pal(n_colors, "Set1")[seq_len(n_sections)]
    color_palette <- stats::setNames(color_values, present_sections)
  } else {
    # Ensure all present sections have a color
    missing_sections <- setdiff(present_sections, names(color_palette))
    if (length(missing_sections) > 0) {
      warning("Some sections missing from color_palette: ", 
              paste(missing_sections, collapse = ", "),
              ".Using default colors for these.")
      n_missing <- length(missing_sections)
      extra_colors <- RColorBrewer::brewer.pal(max(3, n_missing), "Set2")[seq_len(n_missing)]
      color_palette <- c(color_palette, stats::setNames(extra_colors, missing_sections))
    }
  }
  
  # ---------------------------
  # 10. Calculate axis limits
  # ---------------------------
  calc_xlim <- if (! is.null(xlim)) {
    xlim
  } else {
    range(draws_filtered$x, na.rm = TRUE)
  }
  
  breaks <- x_breaks %||% ggplot2::waiver()
  
  # ---------------------------
  # 11. Set up y-axis labels based on prob_reference
  # ---------------------------
  if (prob_reference == "null") {
    y_left_label <- paste0("Probability ", measure, " < ", null_value)
    y_right_label <- paste0("Probability ", measure, " > ", null_value)
  } else {
    # null_range reference
    y_left_label <- paste0("Probability ", measure, " < ", null_range[1])
    y_right_label <- paste0("Probability ", measure, " > ", null_range[2])
  }
  
  # ---------------------------
  # 12. Create ECDF plot
  # ---------------------------
  ecdf_plot <- ggplot2::ggplot(draws_filtered, ggplot2::aes(x = x, color = section_label)) +
    # Add null range shading if requested
    {if (isTRUE(add_null_range) && !is.null(null_range)) {
      ggplot2::annotate("rect",
                        xmin = null_range[1], xmax = null_range[2],
                        ymin = -Inf, ymax = Inf,
                        fill = scales::alpha(color_null_range, 0.2),
                        color = NA)
    }} +
    # Add null range boundary lines
    {if (isTRUE(add_null_range) && !is.null(null_range)) {
      ggplot2::geom_vline(xintercept = null_range,
                          linetype = "dotted", color = "grey30")
    }} +
    # Add ECDF
    ggplot2::stat_ecdf(geom = "step", linewidth = 0.9) +
    # Null line
    ggplot2::geom_vline(xintercept = null_value, linewidth = 0.8, color = "black") +
    # Y-axis scales
    ggplot2::scale_y_continuous(
      name = y_left_label,
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1),
      labels = scales::percent_format(),
      sec.axis = ggplot2::sec_axis(
        ~ 1 - .,
        name = y_right_label,
        breaks = seq(0, 1, 0.1),
        labels = scales::percent_format()
      )
    ) +
    # Color scale with specified order - use limits to control legend
    ggplot2::scale_color_manual(
      values = color_palette,
      limits = present_sections,
      name = "Analysis"
    ) +
    # Theme
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = c(0.01, 0.99),
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_rect(
        fill = "white", color = "grey60", linewidth = 0.3
      ),
      legend.title = ggplot2::element_text(family = font, face = "bold"),
      legend.text = ggplot2::element_text(family = font),
      axis.title = ggplot2::element_text(family = font),
      axis.text = ggplot2::element_text(family = font),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(family = font, face = "bold"),
      plot.subtitle = ggplot2::element_text(family = font),
      # Match x-axis line and ticks to null line style
      axis.line.x.bottom = ggplot2::element_line(linewidth = 0.8, color = "black"),
      axis.ticks.x = ggplot2::element_line(linewidth = 0.8, color = "black"),
      axis.ticks.length.x = ggplot2::unit(0.15, "cm")
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(linetype = "solid", shape = NA)
      )
    )
  
  # Apply appropriate x-axis scale
  if (isTRUE(props$log_scale)) {
    ecdf_plot <- ecdf_plot +
      ggplot2::scale_x_log10(
        breaks = breaks,
        expand = c(0, 0)
      ) +
      ggplot2:: coord_cartesian(xlim = calc_xlim, clip = "off") +
      ggplot2:: labs(x = props$x_label)
  } else {
    ecdf_plot <- ecdf_plot +
      ggplot2::scale_x_continuous(
        breaks = breaks,
        expand = c(0, 0)
      ) +
      ggplot2::coord_cartesian(xlim = calc_xlim, clip = "off") +
      ggplot2::labs(x = props$x_label)
  }
  
  # ---------------------------
  # 13. Create density plot if requested
  # ---------------------------
  if (isTRUE(show_density)) {
    # Modify ECDF to remove x-axis elements
    ecdf_plot <- ecdf_plot +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2:: element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.line.x.bottom = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = 2, l = 5.5, unit = "pt")
      )
    
    # Create density plot - NO FILL, only colored outlines
    density_plot <- ggplot2::ggplot(draws_filtered, ggplot2::aes(x = x, color = section_label)) +
      # Add null range shading if requested
      {if (isTRUE(add_null_range) && !is.null(null_range)) {
        ggplot2::annotate("rect",
                          xmin = null_range[1], xmax = null_range[2],
                          ymin = -Inf, ymax = Inf,
                          fill = scales::alpha(color_null_range, 0.2),
                          color = NA)
      }} +
      # Add null range boundary lines
      {if (isTRUE(add_null_range) && !is.null(null_range)) {
        ggplot2::geom_vline(xintercept = null_range,
                            linetype = "dotted", color = "grey30")
      }} +
      # Add density curves - NO FILL
      ggdist::stat_slab(
        fill = NA,
        slab_linewidth = 0.9,
        normalize = "all"
      ) +
      # Null line
      ggplot2::geom_vline(xintercept = null_value, linewidth = 0.8, color = "black") +
      # Color scale with specified order
      ggplot2::scale_color_manual(
        values = color_palette,
        limits = present_sections,
        name = "Analysis"
      ) +
      # Y-axis with some expansion for labels
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
      # Theme
      ggplot2::theme_light() +
      ggplot2::theme(
        legend.position = "none",
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2:: element_blank(),
        axis.title.x = ggplot2:: element_text(family = font),
        axis.text.x = ggplot2::element_text(family = font),
        panel.grid.minor = ggplot2::element_blank(),
        #panel.border = ggplot2::element_blank(),
        # Match x-axis line and ticks to null line style
        axis.line.x.bottom = ggplot2::element_line(linewidth = 0.8, color = "black"),
        axis.ticks.x = ggplot2::element_line(linewidth = 0.8, color = "black"),
        axis.ticks.length.x = ggplot2::unit(0.15, "cm"),
        plot.margin = ggplot2::margin(t = 2, r = 5.5, b = 5.5, l = 5.5, unit = "pt")
      ) +
      ggplot2::labs(x = props$x_label)
    
    # Apply appropriate x-axis scale to density plot
    if (isTRUE(props$log_scale)) {
      density_plot <- density_plot +
        ggplot2::scale_x_log10(
          breaks = breaks,
          expand = c(0, 0)
        ) +
        ggplot2::coord_cartesian(xlim = calc_xlim, clip = "off")
    } else {
      density_plot <- density_plot +
        ggplot2::scale_x_continuous(
          breaks = breaks,
          expand = c(0, 0)
        ) +
        ggplot2:: coord_cartesian(xlim = calc_xlim, clip = "off")
    }
    
    # Add Favours labels using annotation_custom for proper positioning
    density_plot <- density_plot +
      ggplot2::annotation_custom(
        grid::textGrob(
          label = paste0("Favours\n", label_intervention),
          x = grid::unit(0.05, "npc"),
          y = grid::unit(0.6, "npc"),
          just = c("left", "bottom"),
          gp = grid:: gpar(col = "grey30", fontsize = 9, fontface = "bold.italic", fontfamily = font)
        ), xmin = calc_xlim[1] - 0.01, xmax = calc_xlim[2], ymin = -Inf, ymax = Inf
      ) +
      ggplot2:: annotation_custom(
        grid::textGrob(
          label = paste0("Favours\n", label_control),
          x = grid::unit(0.97, "npc"),
          y = grid::unit(0.6, "npc"),
          just = c("right", "bottom"),
          gp = grid::gpar(col = "grey30", fontsize = 9, fontface = "bold.italic", fontfamily = font)
        ), xmin = calc_xlim[1], xmax = calc_xlim[2]- 0.01, ymin = -Inf, ymax = Inf
      )
    
    # Combine plots using patchwork with axis alignment
    final_plot <- (ecdf_plot /patchwork::plot_spacer()/ density_plot) +
      patchwork::plot_layout(
        ncol = 1,
        heights = c(1, 0.0, 0.35),
        axes = "collect_x",
        axis_titles = "collect_x"
      ) &
      ggplot2::coord_cartesian(xlim = calc_xlim)
    
    return(final_plot)
  }
  
  ecdf_plot
}
