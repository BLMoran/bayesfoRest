#' Internal function to create table for left side of forest plot
#'
#' @noRd
forest.table.left_fn <- function(forest.data.summary,
                                 subgroup = FALSE,
                                 label_control = "Control",
                                 label_intervention = "Intervention",
                                 measure = "OR",
                                 font = NULL) {

  is_continuous <- measure %in% c("MD", "SMD")

  # Define column labels based on effect type
  if (is_continuous) {
    control_label <- paste0(label_control, "\nMean (SD)")
    int_label <- paste0(label_intervention, "\nMean (SD)")
  } else {
    control_label <- paste(label_control, "\n(Events/Total)")
    int_label <- paste(label_intervention, "\n(Events/Total)")
  }

  # Choose correct columns for data
  if (is_continuous) {
    selected_cols <- c("Author", "Year", "N_int", "int_mean_sd", "N_ctrl", "ctrl_mean_sd")
  } else {
    selected_cols <- c("Author", "Year", "int_outcome_frac", "control_outcome_frac")
  }
  if (isTRUE(subgroup)) {
    selected_cols <- c("Author", "Year", "Subgroup", selected_cols[!(selected_cols %in% c("Author", "Year"))])
  }

  df <- forest.data.summary |> dplyr::select(dplyr::any_of(selected_cols))

  df <- df |> dplyr::mutate(
    Author = dplyr::if_else(
      is.na(Year),
      Author,
      paste0(Author, " (", Year, ")")
    )
  )

  # Replace NA values with blanks for 'Pooled Effect' and 'Prediction' rows
  if (is_continuous) {
    df <- df |> dplyr::mutate(
      int_mean_sd  = dplyr::if_else(Author %in% c("Pooled Effect", "Prediction"), NA_character_, int_mean_sd),
      ctrl_mean_sd = dplyr::if_else(Author %in% c("Pooled Effect", "Prediction"), NA_character_, ctrl_mean_sd),
      N_int  = dplyr::if_else(Author == "Prediction" | N_int  == "NA", NA_character_, N_int),
      N_ctrl = dplyr::if_else(Author == "Prediction" | N_ctrl == "NA", NA_character_, N_ctrl)
    )
  }

  # Base gt table
  if (isFALSE(subgroup)) {
    forest.table.left <- df |>
      gt::gt() |>
      gt::cols_label(Author = "Study") |>
      gt::cols_align(align = "left") |>
      gt::tab_style(
        style = gt::cell_text(style = "italic", weight = "bold"),
        locations = gt::cells_body(rows = Author %in% c("Pooled Effect", "Prediction")))
  } else {
    forest.table.left <- df |>
      gt::gt(groupname_col = "Subgroup") |>
      gt::tab_options(row_group.font.weight = "bold") |>
      gt::cols_label(Author = gt::md("Subgroup/ \nStudy")) |>
      gt::cols_align(align = "left") |>
      gt::tab_style(
        gt::cell_text(color = "white"),
        locations = gt::cells_row_groups(groups = "Overall")) |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold", style = "italic"),
        locations = gt::cells_body(rows = Author == "Overall Effect")) |>
      gt::tab_style(
        style = gt::cell_text(style = "italic", weight = "bold", color = "grey60"),
        locations = gt::cells_body(rows = Author == "Pooled Effect")) |>
      # Subgroup-level Prediction: same style as Pooled Effect
      gt::tab_style(
        style = gt::cell_text(style = "italic", weight = "bold", color = "grey60"),
        locations = gt::cells_body(rows = Author == "Prediction" & Subgroup != "Overall")) |>
      # Overall Prediction: same style as Overall Effect
      gt::tab_style(
        style = gt::cell_text(style = "italic", weight = "bold"),
        locations = gt::cells_body(rows = Author == "Prediction" & Subgroup == "Overall")) |>
      gt::tab_style(
        style = gt::cell_text(indent = gt::px(10)),
        locations = gt::cells_body(rows = !Author %in% c("Overall Effect") &
                                     !(Author == "Prediction" & Subgroup == "Overall"),
                                   columns = Author))
  }

  # Column label mapping
  col_labels <- if (is_continuous) {
    c(N_int = "N", int_mean_sd = gt::md(int_label),
      N_ctrl = "N", ctrl_mean_sd = gt::md(control_label))
  } else {
    c(control_outcome_frac = gt::md(control_label),
      int_outcome_frac = gt::md(int_label))
  }

  # Final styling
  forest.table.left <- forest.table.left |>
    gt::cols_label(!!!col_labels) |>
    gt::tab_options(
      column_labels.font.weight = "bold",
      table.border.top.color = "white",
      table.border.bottom.color = "white",
      table.font.names = font) |>
    gt::cols_hide(columns = Year) |>
    gt::opt_table_lines(extent = "none") |>
    gt::tab_style(
      style = gt::cell_fill(color = "grey95"),
      locations = gt::cells_body(rows = Author %in% c("Pooled Effect", "Overall Effect", "Prediction"))
    )

  if (measure %in% c("MD", "SMD")) {
    forest.table.left <- forest.table.left |>
      gt::sub_missing(columns = c(int_mean_sd, ctrl_mean_sd, N_int, N_ctrl), missing_text = "")
  }

  return(forest.table.left)
}


#' Internal function to create table on right of forest plot (including RoB)
#'
#' @param rob_cols Risk of bias columns to show (from `rob_columns()`).
#' @param spec A `rob_spec` object (only used when `add_rob = TRUE`).
#'
#' @noRd
forest.table.right_fn <- function(df,
                                  subgroup = FALSE,
                                  measure = "OR",
                                  add_rob = FALSE,
                                  rob_cols = character(),
                                  spec = NULL,
                                  effect_label = NULL,
                                  font = NULL) {
  if (isFALSE(add_rob)) rob_cols <- character()

  if (isFALSE(subgroup)) {
    forest.table.right <- df |>
      dplyr::select(Author, weighted_effect, unweighted_effect, dplyr::all_of(rob_cols)) |>
      gt::gt() |>
      gt::tab_style(
        style = list(
          gt::cell_text(
            style = "italic",
            weight = "bold")),
        locations = gt::cells_body(
          rows = Author %in% c("Pooled Effect", "Prediction")))

  } else if (isTRUE(subgroup)) {
    forest.table.right <- df |>
      dplyr::select(Author, Subgroup, weighted_effect, unweighted_effect, dplyr::all_of(rob_cols)) |>
      gt::gt(groupname_col = "Subgroup") |>
      gt::tab_style(
        gt::cell_text(color = "white"),
        locations = gt::cells_row_groups(groups = gt::everything())) |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold", style = "italic"),
        locations = gt::cells_body(rows = Author == "Overall Effect")) |>
      gt::tab_style(
        style = list(
          gt::cell_text(
            style = "italic",
            weight = "bold",
            color = "grey60")),
        locations = gt::cells_body(
          rows = Author == "Pooled Effect")) |>
      # Subgroup-level Prediction: grey like Pooled Effect
      gt::tab_style(
        style = list(
          gt::cell_text(
            style = "italic",
            weight = "bold",
            color = "grey60")),
        locations = gt::cells_body(
          rows = Author == "Prediction" & Subgroup != "Overall")) |>
      # Overall Prediction: bold italic, no grey
      gt::tab_style(
        style = list(
          gt::cell_text(
            style = "italic",
            weight = "bold")),
        locations = gt::cells_body(
          rows = Author == "Prediction" & Subgroup == "Overall"))
  }

  forest.table.right <- forest.table.right |>
    gt::tab_options(column_labels.font.weight = "bold") |>
    gt::cols_align(align = "right") |>
    gt::cols_label(.list = create_rob_labels(rob_cols, measure, effect_label)) |>
    gt::tab_options(
      table.border.top.color = "white",
      table.border.bottom.color = "white",
      table.font.names = font) |>
    gt::opt_table_lines(extent = "none") |>
    gt::tab_style(
      style = gt::cell_fill(color = "grey95"),
      locations = gt::cells_body(
        rows = Author %in% c("Pooled Effect", "Overall Effect", "Prediction"),
        columns = c(weighted_effect, unweighted_effect)))

  if (length(rob_cols) > 0) {
    text_cols <- get_rob_text_columns(rob_cols)
    forest.table.right <- forest.table.right |>
      gt::sub_missing(columns = dplyr::all_of(rob_cols), missing_text = "") |>
      gt::cols_align(
        align = "center",
        columns = dplyr::all_of(rob_cols)) |>
      # Reduce horizontal padding only for columns with "Risk", "of", "Bias" text
      gt::tab_style(
        style = "padding-left:0px; padding-right:0px;",
        locations = gt::cells_body(columns = dplyr::all_of(text_cols))) |>
      gt::tab_style(
        style = "padding-left:0px; padding-right:0px;",
        locations = gt::cells_column_labels(columns = dplyr::all_of(text_cols))) |>
      gt::text_transform(
        locations = gt::cells_body(columns = dplyr::all_of(rob_cols)),
        fn = \(x) rob_icon(x, levels = names(spec@levels), height = "1.1em")) |>
      gt::tab_style(
        style = "padding-top:2px; padding-bottom:2px;",
        locations = gt::cells_body(columns = dplyr::all_of(rob_cols))) |>
      gt::tab_style(
        style = "padding-right:2px;",
        locations = gt::cells_body(
          columns = unweighted_effect,
          rows = Author %in% c("Pooled Effect", "Prediction")))
  }

  forest.table.right |>
    gt::cols_hide(columns = Author)
}
