#' Create Density Plot for Sensitivity Analysis
#'
#' @description
#' Creates the central density plot showing posterior distributions
#' across different sensitivity analyses.
#'
#' @param df Data frame containing posterior draws
#' @param measure Effect measure type
#' @param split_color_by_null Logical. Color based on null hypothesis
#' @param color_overall_posterior Fill color for densities
#' @param color_overall_posterior_outline Outline color for densities
#' @param color_favours_control Color when favoring control
#' @param color_favours_intervention Color when favoring intervention
#' @param label_control Label for control group
#' @param label_intervention Label for intervention group
#' @param xlim X-axis limits
#' @param x_breaks X-axis break points
#' @param null_value Null hypothesis value
#' @param null_range Range for practical equivalence
#' @param add_null_range Whether to show null range
#' @param color_null_range Color for null range region
#' @param font Font family
#'
#' @return A ggplot object
#'
#' @keywords internal
#' @noRd
sensitivity_density_plot_fn <- function(df,
                                        measure,
                                        split_color_by_null = FALSE,
                                        color_overall_posterior = "dodgerblue",
                                        color_overall_posterior_outline = "blue",
                                        color_favours_control = "firebrick",
                                        color_favours_intervention = "dodgerblue",
                                        label_control = "Control",
                                        label_intervention = "Intervention",
                                        xlim = NULL,
                                        x_breaks = NULL,
                                        null_value = NULL,
                                        null_range = NULL,
                                        add_null_range = NULL,
                                        color_null_range = "#77bb41",
                                        font = NULL) {
  
  props <- get_measure_properties(measure)
  null_value <- null_value %||% props$null_value
  breaks <- x_breaks %||% ggplot2::waiver()
  
  section_levels <- unique(df$section_label)
  
  df <- df |>
    dplyr::mutate(
      section_label = factor(section_label, levels = section_levels),
      plot_row = factor(
        paste(section_label, "-", prior_label),
        levels = unique(paste(section_label, "-", prior_label))
      ),
      plot_row = forcats::fct_rev(plot_row)
    )
  
  
  # Get the numeric positions of all factor levels
  all_levels <- levels(df$plot_row)
  level_positions <- seq_along(all_levels)
  
  # Find positions where section changes occur
  section_info <- df |>
    dplyr::distinct(plot_row, section_label) |>
    dplyr::arrange(plot_row) |>
    dplyr::mutate(
      level_num = match(plot_row, all_levels),
      section_changes = section_label != dplyr::lag(section_label, default = first(section_label))
    )
  
  # Get positions for lines between sections
  # We want lines after each section ends (at positions 3.5 and 6.5)
  section_breaks <- section_info |>
    dplyr::group_by(section_label) |>
    dplyr::summarise(max_level = max(level_num), .groups = "drop") |>
    dplyr::filter(max_level < max(level_positions)) |>  # Don't add line after last section
    dplyr::pull(max_level) |>
    purrr::map_dbl(~ . + 0.5)  # Add 0.5 to place line after the section
  
  # Position for top line (above the highest level)
  top_line_pos <- max(level_positions) + 0.5
  
  # Fix xlim calculation
  calc_xlim <- if (!is.null(xlim)) xlim else range(df$x, na.rm = TRUE)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(y = plot_row)) +
    {if (isTRUE(add_null_range)) {
      ggplot2::annotate("rect", 
                        xmin = null_range[1], xmax = null_range[2], 
                        ymin = -Inf, ymax = Inf, 
                        fill = scales::alpha(color_null_range, 0.3), 
                        color = NA)
    }} +
    {if (isTRUE(split_color_by_null)) {
      ggdist::stat_slab(
        ggplot2::aes(x = x, fill = ggplot2::after_stat(x > null_value)),
        normalize = "panels",
        height = 0.9,
        colour = color_overall_posterior_outline,
        linewidth = 0.5,  # Thinner outline
        alpha = 0.7,  # More transparent fill
        position = ggplot2::position_nudge(y = -0.5)
      )
    } else {
      ggdist::stat_slab(
        ggplot2::aes(x = x),
        normalize = "panels",
        height = 0.9,
        fill = color_overall_posterior,
        colour = color_overall_posterior_outline,
        linewidth = 0.5,  # Thinner outline
        alpha = 0.7,  # More transparent fill
        position = ggplot2::position_nudge(y = -0.5)
      )
    }} +
    {if (isTRUE(split_color_by_null)) {
      ggplot2::scale_fill_manual(
        values = c(
          "FALSE" = color_favours_intervention,
          "TRUE"  = color_favours_control
        ),
        guide = "none"
      )
    }} +
    
    # Null line
    ggplot2::geom_vline(
      xintercept = null_value,
      linewidth = 1,
      color = "black"
    ) +
    
    # Add horizontal lines between sections (thinner and lighter)
    {if (length(section_breaks) > 0) {
      purrr::map(section_breaks, ~ {
        ggplot2::geom_hline(
          yintercept = .x,
          color = "grey85",  # Very light grey
          linewidth = 0.3    # Thin line
        )
      })
    }} +
    
    # Add grey line at the top of the plot
    ggplot2::geom_hline(
      yintercept = top_line_pos,
      color = "grey60",  # Medium grey for top line
      linewidth = 0.5
    ) +
    
    # Add grey line at the bottom of the plot
    ggplot2::geom_hline(
      yintercept = 0.5,
      color = "grey60",  # Medium grey for bottom line
      linewidth = 0.5
    ) +
    
    # Extend y-axis slightly to show top and bottom lines
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    
    # Add favours labels
    ggplot2::annotation_custom(
      grid::textGrob(
        label = paste(" Favours\n", label_control),  
        x = grid::unit(0.97, "npc"), 
        y = grid::unit(1.02, "npc"), 
        just = c("right", "bottom"),
        gp = grid::gpar(col = "grey30", fontsize = 9, fontfamily = font)
      ),
      xmin = calc_xlim[1], xmax = calc_xlim[2]- 0.01, ymin = -Inf, ymax = Inf
    ) +
    
    ggplot2::annotation_custom(
      grid::textGrob(
        label = paste(" Favours\n", label_intervention),
        x = grid::unit(0.05, "npc"), 
        y = grid::unit(1.02, "npc"), 
        just = c("left", "bottom"),
        gp = grid::gpar(col = "grey30", fontsize = 9, fontfamily = font)
      ),
      xmin = calc_xlim[1] - 0.01, xmax = calc_xlim[2], ymin = -Inf, ymax = Inf
    ) +
    
    ggplot2::coord_cartesian(
      xlim = calc_xlim,
      ylim = c(0.5, max(level_positions) + 0.5),
      clip = "off"
    ) +
    
    
    ggplot2::theme_light() +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(vjust = -0.5, family = font),
      panel.grid   = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(0, 0, 0, 0),
      axis.line.x.bottom = ggplot2::element_line(
        linewidth = 0.5, color = "grey60"),
      axis.line.x.top = ggplot2::element_line(
        linewidth = 0.5, color = "grey60"),
    ) +
    
    ggplot2::labs(x = props$x_label)
  
  # Apply appropriate scale
  if (isTRUE(props$log_scale)) {
    p <- p + ggplot2::scale_x_log10(
      breaks = breaks,
      limits = calc_xlim,
      expand = c(0, 0)
    )
  } else {
    p <- p + ggplot2::scale_x_continuous(
      breaks = breaks,
      limits = calc_xlim,
      expand = c(0, 0)
    )
  }
  
  p
}