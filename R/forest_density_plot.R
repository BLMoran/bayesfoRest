#' Internal function to generate plot section of forest plot
#'
#' @noRd
study.density.plot_fn  <- function(df,
                                   model,
                                   measure = measure,
                                   subgroup = FALSE,
                                   subgroup_order = NULL,
                                   color_palette = NULL,
                                   color_study_posterior = "dodgerblue",
                                   color_study_posterior_outline = "blue",
                                   color_pooled_posterior = "blue",
                                   color_overall_posterior = "darkblue",
                                   color_shrinkage_outline = "purple",
                                   color_pointinterval = "purple",
                                   color_shrinkage_fill = NA,
                                   split_color_by_null = FALSE,
                                   color_favours_control = "firebrick",
                                   color_favours_intervention = "dodgerblue",
                                   label_control = "Control",
                                   label_intervention = "Intervention",
                                   shrinkage_output = "density",
                                   xlim = NULL,
                                   x_breaks = NULL,
                                   null_value = NULL,
                                   add_rope = FALSE,
                                   rope_value = NULL,
                                   rope_color = "grey50",
                                   font = NULL,
                                   reverse_arms = FALSE){
  
  # Filter out "No Pooled Effect" rows at the beginning
  df <- df |> dplyr::filter(Author != "No Pooled Effect")
  
  # Apply color palette if provided
  if (!is.null(color_palette)) {
    # Extract colors from the specified palette
    # Assuming we need at least 6 colors for the various elements
    palette_colors <- paletteer::paletteer_d(color_palette, n = 6)
    
    # Override default colors with palette colors
    color_study_posterior <- palette_colors[2]
    color_study_posterior_outline <- palette_colors[1]
    color_pooled_posterior <- palette_colors[3]
    color_overall_posterior <- palette_colors[5]
    color_shrinkage_outline <- palette_colors[6]
    color_pointinterval <- palette_colors[6]
    
    # If you want to use a fill color from the palette instead of NA
    if (!is.na(color_shrinkage_fill)) {
      color_shrinkage_fill <- scales::alpha(palette_colors[6], 0.3)
    }
  }
  
  # Get effect size properties
  props <- get_measure_properties(measure) 
  
  # Determine x_breaks to use: custom if provided, otherwise use measure defaults
  breaks <- if (!is.null(x_breaks)) x_breaks else ggplot2::waiver()
  
  # Determine null line to use: custom if provided, otherwise use measure default
  null_value <- if (!is.null(null_value)) null_value else props$null_value
  
  # Process ROPE values if provided
  if (isTRUE(add_rope)) {
    if (is.null(rope_value)) {
      # Use Kruschke's recommended defaults if no rope_value provided
      if (measure %in% c("OR", "HR", "RR", "IRR")) {
        # For ratio measures: ROPE around 1.0 with +/-10% range (0.9 to 1.1)
        rope_range <- c(0.9, 1.1)
      } else if (measure == "SMD") {
        # For standardized mean difference: Cohen's small effect size +/-0.1
        rope_range <- c(-0.1, 0.1)
      } else if (measure == "MD") {
        # For mean difference: no universal default, require user specification
        stop("For mean differences, rope_value must be specified as it depends on the measurement scale")
      }
    } else {
      if (length(rope_value) == 1) {
        # Single value: create symmetric range around null line
        rope_range <- c(null_line_to_use - rope_value, null_line_to_use + rope_value)
      } else if (length(rope_value) == 2) {
        # Two values: use as explicit range
        rope_range <- sort(rope_value)  # Ensure correct order
      } else {
        stop("rope_value must be either a single value or a vector of length 2")
      }
    }
  }
  
  # Optimise data structure: separate study-level data from posterior draws
  if (isTRUE(subgroup)){
    study.effects <- df |> dplyr::distinct(Author_ordered, yi, vi) |>
      dplyr::mutate(
        Author_ordered = factor(Author_ordered, levels = as.character(1:max(Author_ordered, na.rm = TRUE))),
        Author = forcats::fct_rev(Author_ordered))
    
    posterior.draws <- df |>
      dplyr::mutate(
        Author_pooled = Author,
        Author_ordered = factor(Author_ordered, levels = as.character(1:max(Author_ordered, na.rm = TRUE))),
        Author = forcats::fct_rev(Author_ordered))
    
  } else if (isFALSE(subgroup)){
    study.effects <- df |> dplyr::distinct(Author, yi, vi)
    posterior.draws <- df
  }
  
  if (isTRUE(props$log_scale)){
    x.min <- 0.1
    x.max <- ceiling(max(exp(df$yi + 1.96 * sqrt(df$vi)), na.rm = TRUE))
  } else if (isFALSE(props$log_scale)){
    x.min <- floor(min(df$yi - 1.96 * sqrt(df$vi), na.rm = TRUE))
    x.max <- ceiling(max(df$yi + 1.96 * sqrt(df$vi), na.rm = TRUE))
  }
  
  # Set xlim to given value or calculated range
  if (!is.null(xlim)){
    calc_xlim <- xlim
  } else {
    calc_xlim <- c(x.min, x.max)
  }
  
  if (isTRUE(props$log_scale)){
    posterior.draws <- posterior.draws |> 
      dplyr::mutate(
        x_studies = exp(b_Intercept))
    
    study.effects <- study.effects |> 
      dplyr::mutate(
        xdist = exp(distributional::dist_normal(mean = yi, sd = sqrt(vi))))
  } else {
    posterior.draws <- posterior.draws |> 
      dplyr::mutate(
        x_studies = b_Intercept)
    
    study.effects <- study.effects |> 
      dplyr::mutate(
        xdist = distributional::dist_normal(mean = yi, sd = sqrt(vi)))
  }
  
  study.density.plot <-
    ggplot2::ggplot(ggplot2::aes(y = Author), data = posterior.draws) +
    # Add ROPE shading first (so it appears behind other elements)
    {if (isTRUE(add_rope)) {
      ggplot2::annotate("rect", 
                        xmin = rope_range[1], xmax = rope_range[2], 
                        ymin = -Inf, ymax = Inf, 
                        fill = scales::alpha(rope_color, 0.3), 
                        color = NA)
    }} +
    {if (isTRUE(split_color_by_null)) {
      list(
        ggdist::stat_slab(
          ggplot2::aes(xdist = xdist, fill = ggplot2::after_stat(x > null_value)),
          slab_linewidth = 0.5, alpha = 0.7, limits = calc_xlim, height = 0.9,
          data = study.effects, colour = color_study_posterior_outline))
    } else {
      ggdist::stat_slab(
        ggplot2::aes(xdist = xdist),
        slab_linewidth = 0.5, alpha = 0.7, limits = calc_xlim, height = 0.9,
        data = study.effects, colour = color_study_posterior_outline,
        fill = color_study_posterior)
    }} +
    {if (isTRUE(split_color_by_null)) {
      list(
        ggdist::stat_slab(
          ggplot2::aes(
            x = x_studies, y = Author, fill = ggplot2::after_stat(x > null_value)),
          data = posterior.draws |>
            dplyr::filter(if (isTRUE(subgroup)) {
              Author_pooled == "Pooled Effect"
            } else {
              Author == "Pooled Effect"
            }),
          height = 0.9, normalize = "panels", colour = color_study_posterior_outline))
    } else {
      ggdist::stat_slab(
        ggplot2::aes(x = x_studies, y = Author),
        data = posterior.draws |>
          dplyr::filter(if (isTRUE(subgroup)) {
            Author_pooled == "Pooled Effect"
          } else {
            Author == "Pooled Effect"
          }),
        fill = color_pooled_posterior, height = 0.9, normalize = "panels")
    }} +
    # Add overall effect slab when subgroup = TRUE
    {if (isTRUE(subgroup)) {
      if (isTRUE(split_color_by_null)) {
        list(
          ggdist::stat_slab(
            ggplot2::aes(
              x = x_studies, y = Author, fill = ggplot2::after_stat(x > null_value)),
            data = posterior.draws |> dplyr::filter(Author_pooled == "Overall Effect"),
            height = 0.9, normalize = "panels", colour = color_study_posterior_outline))
      } else {
        ggdist::stat_slab(
          ggplot2::aes(x = x_studies, y = Author),
          data = posterior.draws |> dplyr::filter(Author_pooled == "Overall Effect"),
          fill = color_overall_posterior, height = 0.9, normalize = "panels")
      }
    }} +
    {if (isTRUE(split_color_by_null)) {
      ggplot2::scale_fill_manual(
        values = c(
          "FALSE" = if (reverse_arms) color_favours_control else color_favours_intervention,
          "TRUE"  = if (reverse_arms) color_favours_intervention else color_favours_control),  
        guide = "none")
    }} +
    # Add null value
    ggplot2::geom_vline(xintercept = null_value, color = "black", linewidth = 1) +
    ggplot2::coord_cartesian(xlim= calc_xlim, clip = "off") +
    ggplot2::theme_light() +
    ggplot2::theme(
      axis.ticks.y = ggplot2::element_blank(), 
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(vjust=-0.5, family = font), 
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0,0,0,0), 
      panel.border = ggplot2::element_blank(),
      axis.line.x.top = ggplot2::element_line(color = "grey60", linewidth = 0.75),
      axis.line.x.bottom = ggplot2::element_line(color = "black", linewidth = 0.75),
      axis.text.x.top = ggplot2::element_blank(), 
      axis.ticks.x.top = ggplot2::element_blank(),
      axis.text.x.bottom = ggplot2::element_text(colour = "black", family = font)) +
    ggplot2::guides(x.sec = "axis", y.sec = "axis") +
    ggplot2::annotation_custom(grid::textGrob(
      label = paste(" Favours\n", if (reverse_arms) label_intervention else label_control),  
      x = grid::unit(1, "npc"), y = grid::unit(1.02, "npc"), just = c("right", "bottom"),
      gp = grid::gpar(col = "grey30", fontsize = 10, fontfamily = font)),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    ggplot2::annotation_custom(grid::textGrob(
      label = paste(" Favours\n", if (reverse_arms) label_control else label_intervention),
      x = grid::unit(0, "npc"), y = grid::unit(1.02, "npc"), just = c("left", "bottom"),
      gp = grid::gpar(col = "grey30", fontsize = 10, fontfamily = font)),
      xmin = calc_xlim[1] - 0.01, xmax = calc_xlim[2], ymin = -Inf, ymax = Inf) +
    ggplot2::labs(x= props$x_label) +
    ggplot2::ylab(NULL) +
    ggplot2::geom_hline(yintercept = 1, color = "black", linewidth = 0.75)
  
  if (isTRUE(props$log_scale)){
    if (isTRUE(reverse_arms)) {
      study.density.plot <- study.density.plot +
        ggplot2::scale_x_continuous(
          transform = scales::compose_trans("log10", "reverse"),
          breaks = breaks, expand = c(0, 0), limits = calc_xlim) +
        ggplot2::geom_vline(xintercept = exp(brms::fixef(model)[1, 1]), color = "grey60", linewidth = 1) +
        ggplot2::geom_vline(xintercept = exp(brms::fixef(model)[1, 3:4]), color = "grey60", linetype = 2)
    } else {
      study.density.plot <- study.density.plot +
        ggplot2::scale_x_log10(breaks = breaks, expand = c(0, 0), limits = calc_xlim)+
        ggplot2::geom_vline(xintercept = exp(brms::fixef(model)[1, 1]), color = "grey60", linewidth = 1) +
        ggplot2::geom_vline(xintercept = exp(brms::fixef(model)[1, 3:4]), color = "grey60", linetype = 2)
    }
  } else {
    if (isTRUE(reverse_arms)) {
      study.density.plot <- study.density.plot +
        ggplot2::scale_x_continuous(
          transform = scales::transform_reverse(),
          breaks = breaks, expand = c(0, 0), limits = calc_xlim) +
        ggplot2::geom_vline(xintercept = brms::fixef(model)[1, 1], color = "grey60", linewidth = 1) +
        ggplot2::geom_vline(xintercept = brms::fixef(model)[1, 3:4], color = "grey60", linetype = 2)
    } else {
      study.density.plot <- study.density.plot +
        ggplot2::scale_x_continuous(breaks = breaks, expand = c(0, 0), limits = calc_xlim) +
        ggplot2::geom_vline(xintercept = brms::fixef(model)[1, 1], color = "grey60", linewidth = 1) +
        ggplot2::geom_vline(xintercept = brms::fixef(model)[1, 3:4], color = "grey60", linetype = 2)
    }
  }
  
  if(isFALSE(subgroup)){
    study.density.plot <- study.density.plot + ggplot2::scale_y_discrete(expand = c(0, 0), limits = rev)
  } else if (isTRUE(subgroup)){
    study.density.plot <- study.density.plot + ggplot2::scale_y_discrete(expand = c(0, 0))
  }
  
  if(is.null(shrinkage_output))"density" else shrinkage_output
  
  if (shrinkage_output == "density") {
    study.density.plot <- study.density.plot +
      ggdist::stat_slab(
        ggplot2::aes(x = x_studies, y = Author),
        data = posterior.draws |> dplyr::filter(if (isTRUE(subgroup)) {Author_pooled != "Pooled Effect" & Author_pooled != "Overall Effect"}
                                                else {Author != "Pooled Effect"}),
        linewidth = 0.5, scale = 0.6, normalize = "panels",
        color = color_shrinkage_outline, fill = color_shrinkage_fill, limits = calc_xlim)
  } else if (shrinkage_output == "pointinterval") {
    study.density.plot <- study.density.plot +
      ggdist::stat_pointinterval(
        ggplot2::aes(x = x_studies, y = Author),
        data = posterior.draws |> dplyr::filter(if (isTRUE(subgroup)) {Author_pooled != "Pooled Effect" & Author_pooled != "Overall Effect"}
                                                else {Author != "Pooled Effect"}),
        linewidth = 1, size = 1, color = color_pointinterval, limits = calc_xlim)
  }
  
  return(study.density.plot)
}