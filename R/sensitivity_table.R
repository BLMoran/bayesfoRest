#' Create Left Table for Sensitivity Plot
#'
#' @description
#' Creates the left-side table showing prior specifications for the
#' sensitivity analysis plot.
#'
#' @param df Data frame containing sensitivity analysis results
#' @param font Optional font family for the table
#' @param math_font Font for mathematical symbols. Default is "STIX Two Math"
#'
#' @return A gt table object
#'
#' @keywords internal
#' @noRd
sensitivity_table_left <- function(
    df,
    font = NULL,
    math_font = "STIX Two Math"
) {
  
  math_cols <- c("mu_prior_unicode", "tau_prior_unicode")
  
  gt_tbl <- df |>
    dplyr::select(section_label, prior_label, mu_prior_unicode, tau_prior_unicode) |>
    dplyr::group_by(section_label) |>
    dplyr::mutate(.last_in_group = dplyr::row_number() == dplyr::n()) |>
    dplyr::ungroup() |> 
    gt::gt(groupname_col = "section_label") |>
    gt::tab_stubhead(label = "Meta-analysis") |>
    gt::cols_label(
      prior_label = "Prior type",
      mu_prior_unicode  = paste0("\u03BC", " Prior"),
      tau_prior_unicode = paste0("\u03C4", " Prior")) |>
    gt::cols_align(align = "left") |>
    gt::tab_options(
      row_group.as_column = TRUE,
      row_group.font.weight = "bold",
      
      column_labels.font.weight = "bold",
      column_labels.border.bottom.color = "grey60",
      column_labels.border.bottom.width = gt::px(2),
      
      table.border.top.color = "white",
      table.border.bottom.color = "grey60",
      table.border.bottom.width = gt::px(2)) |>
    gt::tab_style(
      style = gt::cell_borders(
        sides = "bottom",
        color = "grey60",
        weight = gt::px(2)),
      locations = gt::cells_body(rows = .last_in_group))
  
  
  if (!is.null(font)) {
    gt_tbl <- gt_tbl |>
      gt::tab_style(
        style = gt::cell_text(font = font),
        locations = list(
          gt::cells_body(columns = -dplyr::all_of(math_cols)),
          gt::cells_column_labels(columns = -dplyr::all_of(math_cols)),
          gt::cells_stub(),
          gt::cells_row_groups()))
  }
  
  gt_tbl |>
    gt::tab_style(
      style = gt::cell_text(font = math_font),
      locations = gt::cells_body(columns = dplyr::all_of(math_cols))) |>
    gt::cols_hide(columns = ".last_in_group")
}

#' Create Right Table for Sensitivity Plot
#'
#' @description
#' Creates the right-side table showing effect estimates and probabilities
#' for the sensitivity analysis plot.
#'
#' @param df Data frame containing sensitivity analysis results
#' @param measure Effect measure type
#' @param add_probs Logical. Whether to include probability columns
#' @param font Optional font family for the table
#'
#' @return A gt table object
#'
#' @keywords internal
#' @noRd
sensitivity_table_right <- function(
    df,
    measure,
    add_probs = FALSE,
    font = NULL
) {
  
  pr_cols <- c(
    "pr_benefit",
    "pr_no_benefit",
    "pr_benefit_null_range",
    "pr_no_benefit_null_range"
  )
  
  df <- df |>
    dplyr::mutate(
      estimate = sprintf("%.3f  [%.3f, %.3f]", median, l95, u95)) |>
    dplyr::select(section_label, prior_label, estimate, dplyr::all_of(pr_cols)) |> 
    dplyr::group_by(section_label) |>
    dplyr::mutate(.last_in_group = dplyr::row_number() == dplyr::n()) |>
    dplyr::ungroup()
  
  table_right <- df |> 
    gt::gt() |>
    gt::cols_label(
      estimate = paste(measure, "[95% CrI]"),
      pr_benefit = "Pr(Benefit)",
      pr_no_benefit = "Pr(Harm)",
      pr_benefit_null_range = paste("Pr(Benefit>", "\u03B4", ")"),
      pr_no_benefit_null_range = paste("Pr(Harm>", "\u03B4", ")")) |>
    #gt::cols_align(align = "left", columns = "estimate") |>
    gt::cols_align(align = "center", columns = everything()) |>
    gt::tab_options(
      table.font.names = font,
      column_labels.font.weight = "bold",
      column_labels.border.bottom.color = "grey60",
      column_labels.border.bottom.width = gt::px(2),
      
      row_group.font.weight = "bold",
      row_group.font.size = gt::px(13),
      
      table.border.top.color = "white",
      table.border.bottom.color = "grey60",
      table.border.bottom.width = gt::px(2)) |>
    gt::fmt_percent(
      columns = dplyr::all_of(pr_cols),
      decimals = 1,
      scale_values = FALSE) |>
    gt::tab_style(
      style = gt::cell_borders(
        sides = "bottom",
        color = "grey60",
        weight = gt::px(2)),
      locations = gt::cells_body(
        rows = .last_in_group)) |> 
    gt::cols_hide(columns = c("section_label", "prior_label", ".last_in_group"))
  
  if (isFALSE(add_probs)) {
    table_right <- table_right |> 
      gt::cols_hide(columns = c(dplyr::all_of(pr_cols)))
  }
  
  return(table_right)
}