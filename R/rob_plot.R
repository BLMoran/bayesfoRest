rob_plot <- function(data,
                     studyvar = NULL,
                     sort_studies_by = "author",
                     subgroup = FALSE,
                     sort_subgroup_by = "alphabetical",
                     rob_tool = c("rob2", "robins_i", "robins_e", "quadas2"),
                     add_rob_legend = FALSE,
                     font = NULL,
                     title = NULL,
                     title_align = "left",
                     subtitle = NULL) {
  
  # Input validation
  if (is.null(data)) {
    stop("Data must be provided")
  }
  
  # Determine RoB variables in data
  rob_vars <- c("D1", "D2", "D3", "Overall")
  missing_vars <- setdiff(rob_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop("Risk of Bias columns must be provided for addition to the forest plot")
  }
  
  if (!sort_studies_by %in% c("author", "year", "effect")) {
    stop("sort_studies_by must be one of 'author', 'year', or 'effect'")
  }
  
  # Handle different workflows for single vs subgroup plots
  if (!subgroup) {
    # Single RoB Plot Workflow
    df <- data |> 
      dplyr::select(Author, Year, D1:Overall) |> 
      sort.studies.fn(sort_studies_by)
    
  } else {
    # Subgroup RoB Plot Workflow
    df <- data |> 
      dplyr::select(Author, Year, Subgroup, D1:Overall) |> 
      sort.studies.fn(sort_studies_by)
    
    # Apply subgroup ordering
    if ("Subgroup" %in% names(df)) {
      if (length(sort_subgroup_by) == 1 && sort_subgroup_by == "alphabetical") {
        subgroup_order <- sort(unique(df$Subgroup[!is.na(df$Subgroup)]))
      } else if (is.character(sort_subgroup_by)) {
        # User supplied a custom group order directly
        subgroup_order <- sort_subgroup_by
      } else {
        stop("Invalid input to sort_subgroup_by. Must be 'alphabetical', or a character vector of subgroup names.")
      }
      
      df <- df |>
        dplyr::mutate(Subgroup = factor(Subgroup, levels = subgroup_order)) |>
        dplyr::arrange(Subgroup)
    }
  }
  
  # Create the main RoB table
  rob.table <- rob.table_fn(df, subgroup = subgroup, font = font)
  
  # Add legend if requested
  if (isTRUE(add_rob_legend)) {
    rob_legend <- rob_legend_fn(rob_tool, font)
    rob.table <- patchwork::wrap_table(rob.table, space = "fixed") +
      patchwork::wrap_table(rob_legend, space = "fixed")
  } else {
    rob.table <- patchwork::wrap_table(rob.table, space = "fixed")
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
    
    rob.table <- rob.table + 
      patchwork::plot_annotation(
        title = title,
        subtitle = subtitle,
        theme = title_theme
      )
  }
  
  return(rob.table)
}



rob.table_fn <- function(df,
                              subgroup = FALSE,
                              font = NULL) {
  
  df <- df |> dplyr::mutate(
    Author = dplyr::if_else(
      is.na(Year),
      Author,
      paste0(Author, " (", Year, ")")
    )
  )
  
  # Base gt table
  if (isFALSE(subgroup)) {
    rob.table <- df |>
      gt::gt() |>
      gt::cols_label(Author = "Study")
    
  } else {
    rob.table <- df |>
      gt::gt(groupname_col = "Subgroup") |>
      gt::tab_options(row_group.font.weight = "bold") |>
      gt::tab_style(
        style = gt::cell_fill(color = "grey95"),
        locations = gt::cells_row_groups()) |> 
      gt::cols_label(Author = gt::md("Subgroup/ \nStudy")) |>
      gt::cols_align(align = "left") |>
      gt::tab_style(
        style = gt::cell_text(indent = gt::px(10)),
        locations = gt::cells_body(columns = Author))
  }

  # Final styling
  rob.table <- rob.table |>
    gt::cols_align(align = "left",
                   columns = Author) |> 
    gt::cols_align(align = "center",
                   columns = c(D1:Overall)) |> 
    gt::tab_options(
      column_labels.font.weight = "bold",
      table.border.top.color = "white",
      table.border.bottom.color = "white",
      table.font.names = font) |>
    gt::cols_hide(columns = Year) |>
    gt::opt_table_lines(extent = "none") |> 
    gt::text_case_match(
      "High" ~ fontawesome::fa(name = "circle-plus", fill = "#e32400", height = "1.1em"),
      "Low" ~ fontawesome::fa(name = "circle-minus", fill = "#77bb41", height = "1.1em"),
      "Some concerns" ~ fontawesome::fa(name = "circle-question", fill = "#f5ec00", height = "1.1em"))
  
  return(rob.table)
}
