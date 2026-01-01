#' Internal function to join sections of forest plot
#'
#' @noRd
patchwork_fn <- function(table.left,
                         study.density.plot,
                         table.right,
                         add_rob_legend = FALSE,
                         rob_tool = "rob2",
                         plot_width = NULL,
                         title = NULL,
                         subtitle = NULL,
                         title_align = "center",
                         font = NULL){
  if (is.null(plot_width)){
    plot_width <- 4
  }
  
  # Create the base plot layout
  if (isFALSE(add_rob_legend)) {
    base_plot <- patchwork::wrap_table(table.left, space = "fixed") + 
      study.density.plot + 
      patchwork::wrap_table(table.right, space = "fixed") +
      patchwork::plot_layout(widths = grid::unit(c(-1, plot_width, -1), 
                                                 c("null", "cm", "null")))
  } else {
    rob_legend <- rob_legend_fn(rob_tool, font)
    # Include rob_legend to the right of table.right
    base_plot <- patchwork::wrap_table(table.left, space = "fixed") + 
      study.density.plot + 
      patchwork::wrap_table(table.right, space = "fixed") +
      patchwork::wrap_table(rob_legend, space = "fixed") +
      patchwork::plot_layout(widths = grid::unit(c(-1, plot_width, -1, -1), 
                                                 c("null", "cm", "null", "null")))
  }
  
  # Add title and/or subtitle if provided
  if (!is.null(title) || !is.null(subtitle)) {
    # Set alignment value based on title_align parameter
    hjust_val <- switch(title_align,
                        "left" = 0,
                        "center" = 0.5,
                        "centre" = 0.5,
                        "right" = 1,
                        0.5) # default to center if invalid input
    
    # Create theme elements
    title_theme <- ggplot2::theme()
    
    if (!is.null(title)) {
      title_theme <- title_theme + 
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = 16, 
            face = "bold", 
            hjust = hjust_val,
            margin = ggplot2::margin(b = ifelse(is.null(subtitle), 10, 5)),
            family = font))
    }
    
    if (!is.null(subtitle)) {
      title_theme <- title_theme + 
        ggplot2::theme(
          plot.subtitle = ggplot2::element_text(
            size = 14, 
            hjust = hjust_val,
            margin = ggplot2::margin(b = 10),
            color = "gray30",
            family = font))
    }
    
    base_plot <- base_plot + 
      patchwork::plot_annotation(
        title = title,
        subtitle = subtitle,
        theme = title_theme
      )
  }
  
  return(base_plot)
}