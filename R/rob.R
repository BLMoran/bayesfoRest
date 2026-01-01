#' Internal function to create risk of bias legend table
#'
#' @noRd
rob_legend_fn <- function(rob_tool = c("rob2", "robins_i", "quadas2", "robins_e", "nos"),
                          font = NULL) {
  
  # Get domains for the specified tool
  domains <- get_rob_domains(rob_tool, full_text = TRUE)
  
  # Remove overall domain
  domains <- domains[names(domains) != "Overall"]
  
  # Create the data for the table
  legend_data <- data.frame(
    Code = paste0("(", names(domains), ")"),
    Description = unname(domains),
    stringsAsFactors = FALSE)
  
  # Add empty row for spacing
  legend_data <- rbind(
    legend_data,
    data.frame(Code = "", Description = "", stringsAsFactors = FALSE))
  
  # Add risk level rows with placeholder text
  risk_rows <- data.frame(
    Code = c(" ", "Low", "Some concerns", "High"),
    Description = c(" ", "Low risk", "Some concerns", "High risk"),
    stringsAsFactors = FALSE)
  
  legend_data <- rbind(legend_data, risk_rows)
  
  # Create the gt table
  rob_legend_table <- legend_data |>
    gt::gt() |>
    # Remove column headers
    gt::cols_label(
      Code = "",
      Description = "") |>
    # Add title
    gt::tab_header(
      title = gt::html("<b>Risk of bias legend<b>")) |>
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold", size = gt::px(16))),
      locations = gt::cells_title(groups = "title")) |> 
    # Style the domain codes (first n rows)
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold", size = gt::px(12))),
      locations = gt::cells_body(
        columns = Code,
        rows = 1:length(domains))) |>
    # Transform the risk level text to icons
    gt::text_transform(
      locations = gt::cells_body(
        columns = Code,
        rows = (length(domains) + 2):nrow(legend_data)),
      fn = function(x) {
        dplyr::case_when(
          x == "High" ~ as.character(fontawesome::fa(name = "circle-plus", fill = "#e32400", height = "1.5em")),
          x == "Low" ~ as.character(fontawesome::fa(name = "circle-minus", fill = "#77bb41", height = "1.5em")),
          x == "Some concerns" ~ as.character(fontawesome::fa(name = "circle-question", fill = "#f5ec00", height = "1.5em")),
          TRUE ~ x)
      }) |>
    # Style the risk level descriptions to be bold and italic
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold", style = "italic")),
      locations = gt::cells_body(
        columns = Description,
        rows = (length(domains) + 2):nrow(legend_data))) |>
    # Reduce space in between RoB and legend
    gt::tab_style(
      style = "padding-left:2px;",
      locations = gt::cells_body(
        columns = Code, 
        rows = gt::everything())) |> 
    # Set column widths
    gt::cols_width(
      Code ~ gt::px(40),
      Description ~ gt::px(300)) |>
    # Table options
    gt::tab_options(
      table.border.top.color = "white",
      table.border.bottom.color = "white",
      heading.align = "left",
      heading.padding = gt::px(6),
      table.font.size = gt::px(12),
      table.font.names = font) |>
    gt::opt_table_lines(extent = "none")
  
  return(rob_legend_table)
}


#' Internal function to get RoB domains for RoB legend
#'
#' @noRd
get_rob_domains <- function(rob_tool = c("rob2", "robins_i", "quadas2", "robins_e"),
                            full_text = TRUE) {
  rob_tool <- match.arg(rob_tool)
  
  domains <- switch(rob_tool,
                    rob2 = c(
                      D1 = "Bias arising from the randomization process",
                      D2 = "Bias due to deviations from intended interventions",
                      D3 = "Bias due to missing outcome data",
                      D4 = "Bias in measurement of the outcome",
                      D5 = "Bias in selection of the reported result",
                      Overall = "Overall risk of bias"
                    ),
                    robins_i = c(
                      D1 = "Bias due to confounding",
                      D2 = "Bias in selection of participants into the study",
                      D3 = "Bias in classification of interventions",
                      D4 = "Bias due to deviations from intended interventions",
                      D5 = "Bias due to missing data",
                      D6 = "Bias in measurement of outcomes",
                      D7 = "Bias in selection of the reported result",
                      Overall = "Overall risk of bias"
                    ),
                    quadas2 = c(
                      D1 = "Patient selection",
                      D2 = "Index test(s)",
                      D3 = "Reference standard",
                      D4 = "Flow and timing",
                      Overall = "Overall risk of bias"
                    ),
                    robins_e = c(
                      D1 = "Bias due to confounding",
                      D2 = "Bias in selection of participants into the study",
                      D3 = "Bias in classification of exposures",
                      D4 = "Bias due to departures from intended exposures",
                      D5 = "Bias due to missing data",
                      D6 = "Bias in measurement of outcomes",
                      D7 = "Bias in selection of the reported result",
                      Overall = "Overall risk of bias"
                    )
  )
  
  if (isTRUE(full_text)){ 
    return(domains)
  } else {
    return(names(domains))
  }
}

#' Internal function to get columns that contain Risk/of/Bias text
#'
#' @noRd
get_rob_text_columns <- function(df) {
  
  # Get the column names that start with D or are "Overall"
  rob_cols <- names(df)[grepl("^D\\d+$|^Overall$", names(df))]
  
  # Split "Risk of Bias" across the middle columns
  rob_text <- c("Risk", "of", "Bias")
  n_rob_cols <- length(rob_cols)
  
  # Calculate which columns get the text
  if (n_rob_cols >= 3) {
    start_pos <- ceiling((n_rob_cols - 3) / 2) + 1
    text_positions <- start_pos:(start_pos + 2)
  } else {
    text_positions <- 1:min(n_rob_cols, 3)
  }
  
  # Return only the columns that have the Risk/of/Bias text
  return(rob_cols[text_positions])
}

#' Internal function to create RoB labels based on how many domains are present
#'
#' @noRd
create_rob_labels <- function(df, measure) {
  
  # Get the column names that start with D or are "Overall"
  rob_cols <- names(df)[grepl("^D\\d+$|^Overall$", names(df))]
  
  # Create base labels
  base_labels <- list(
    weighted_effect = gt::md(paste("Shrinkage", measure, "\n", "[95% CrI]")),
    unweighted_effect = gt::md(paste("Observed", measure, "\n", "[95% CrI]"))
  )
  
  # Split "Risk of Bias" across the middle columns
  rob_text <- c("Risk", "of", "Bias")
  n_rob_cols <- length(rob_cols)
  
  # Calculate which columns get the text
  if (n_rob_cols >= 3) {
    start_pos <- ceiling((n_rob_cols - 3) / 2) + 1
    text_positions <- start_pos:(start_pos + 2)
  } else {
    text_positions <- 1:min(n_rob_cols, 3)
    rob_text <- rob_text[1:length(text_positions)]
  }
  
  # Create labels for ROB columns using purrr
  rob_labels <- purrr::imap(rob_cols, ~ {
    col_name <- .x
    position <- .y
    
    if (position %in% text_positions) {
      text_idx <- which(text_positions == position)
      gt::md(paste0("**", rob_text[text_idx], "**\n", col_name))
    } else {
      gt::md(paste0("&nbsp;\n", col_name))
    }
  }) |>
    purrr::set_names(rob_cols)
  
  # Combine base labels with ROB labels
  return(c(base_labels, rob_labels))
}
