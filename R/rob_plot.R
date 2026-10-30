#' Create Risk of Bias Plot
#'
#' Creates a publication-ready risk of bias plot using gt tables, for single
#' studies or subgroup analyses. Uses the same risk of bias internals as the
#' robbR package.
#'
#' @param data A data frame with a study, year and domain columns for the
#'   chosen `rob_tool` (see Details). Judgements are standardised (e.g.
#'   "low risk" -> "Low"); `NA` -> "Not applicable". Judgements not recognised
#'   by the tool are shown as blank cells, with a warning.
#' @param studyvar,yearvar <[`data-masking`][dplyr::dplyr_data_masking]> Study
#'   label and year columns.
#' @param subgroupvar <[`data-masking`][dplyr::dplyr_data_masking]> Subgroup
#'   column, used when `subgroup = TRUE`.
#' @param effect <[`data-masking`][dplyr::dplyr_data_masking]> Effect size
#'   column, used when sorting by `"effect"`.
#' @param rob_tool Risk of bias tool: `"rob2"` (default), `"robins_i"`,
#'   `"robins_e"`, `"quadas2"`, `"clarity"`, `"rob_nma"`, `"probast"`, or
#'   `"probast_ai"`.
#' @param sort_studies_by How to sort studies: `"author"` (default), `"year"`,
#'   or `"effect"`.
#' @param subgroup Logical. Create a subgroup plot? (default: `FALSE`).
#' @param sort_subgroup_by Subgroup ordering: `"alphabetical"` (default),
#'   `"effect"`, or a character vector of subgroup names in the desired order.
#' @param show_overall Logical. Show the overall column(s)? (default: `TRUE`).
#'   Overall columns absent from `data` are derived from the domains (worst judgement).
#' @param derive_overall Logical. Always re-derive overall column(s) from the
#'   domains, overwriting supplied values? (default: `FALSE`).
#' @param add_rob_legend Logical. Add a legend? (default: `FALSE`).
#' @param study_label Header for the study column. Defaults to `"Trial"` for
#'   randomised trial tools (`"rob2"`, `"clarity"`), `"Review"` for
#'   `"rob_nma"`, and `"Study"` otherwise.
#' @param font Font family.
#' @param title,subtitle Plot title and subtitle.
#' @param title_align Title alignment: `"left"` (default), `"center"`/`"centre"`,
#'   or `"right"`.
#'
#' @details
#' Domain columns by tool (optional columns in brackets are shown only if
#' present; `NA` rows are shown as "Not applicable"):
#' \describe{
#'   \item{rob2}{D1, (D1b), D2-D5, (DS), Overall.}
#'   \item{robins_i}{ROBINS-I V2: D1-D6, Overall.}
#'   \item{robins_e}{D1-D7, Overall.}
#'   \item{quadas2}{D1-D4, Overall.}
#'   \item{clarity}{D1-D10, (D11), (D12). No overall judgement.}
#'   \item{rob_nma}{(R1-R3), D1-D3, Results, Conclusions.}
#'   \item{probast, probast_ai}{D1-D4, (A1-A3), Overall, Applicability.}
#' }
#'
#' @return A patchwork composition of the risk of bias table (and legend).
#'
#' @examples
#' \dontrun{
#' rob_data <- data.frame(
#'   Author = c("Smith", "Jones", "Brown"),
#'   Year = c(2020, 2021, 2022),
#'   D1 = c("Low", "High", "Some concerns"),
#'   D2 = c("Low", "Low", "High"),
#'   D3 = c("Some concerns", "Low", "Low"),
#'   D4 = "Low",
#'   D5 = "Low")
#'
#' rob_plot(rob_data, add_rob_legend = TRUE)
#' }
#'
#' @export
rob_plot <- function(data,
                     studyvar = Author,
                     yearvar = Year,
                     subgroupvar = Subgroup,
                     effect = yi,
                     rob_tool = "rob2",
                     sort_studies_by = c("author", "year", "effect"),
                     subgroup = FALSE,
                     sort_subgroup_by = "alphabetical",
                     show_overall = TRUE,
                     derive_overall = FALSE,
                     add_rob_legend = FALSE,
                     study_label = NULL,
                     font = NULL,
                     title = NULL,
                     title_align = c("left", "center", "centre", "right"),
                     subtitle = NULL) {
  sort_studies_by <- rlang::arg_match(sort_studies_by)
  title_align <- rlang::arg_match(title_align)
  spec <- as_rob_spec(rob_tool)
  uses_effect <- sort_studies_by == "effect" || identical(sort_subgroup_by, "effect")
  data <- dplyr::rename(data, Author = {{ studyvar }}, Year = {{ yearvar }})
  if (isTRUE(subgroup)) data <- dplyr::rename(data, Subgroup = {{ subgroupvar }})
  if (uses_effect) data <- dplyr::rename(data, yi = {{ effect }})
  data <- prepare_rob_data(data, spec, derive_overall = derive_overall)
  rob_cols <- rob_columns(data, spec)
  if (isFALSE(show_overall)) rob_cols <- setdiff(rob_cols, names(spec@overall))
  df <- data |>
    dplyr::select("Author", "Year", dplyr::any_of(c(if (isTRUE(subgroup)) "Subgroup", if (uses_effect) "yi")), dplyr::all_of(rob_cols)) |>
    sort_studies_fn(sort_studies_by) |>
    order_subgroups(subgroup, sort_subgroup_by) |>
    dplyr::select(-dplyr::any_of("yi"))
  plot <- rob_table_fn(df, spec, rob_cols, subgroup = subgroup, study_label = study_label, font = font) |>
    patchwork::wrap_table(space = "fixed")
  if (isTRUE(add_rob_legend)) {
    levels_present <- df |>
      dplyr::select(dplyr::all_of(rob_cols)) |>
      purrr::map(as.character) |>
      purrr::list_c() |>
      unique()
    legend <- rob_legend_fn(spec, domains_present = dplyr::intersect(rob_cols, names(spec@domains)), levels_present = levels_present, font = font)
    plot <- plot + patchwork::wrap_table(legend, space = "fixed")
  }
  add_rob_titles(plot, title, subtitle, title_align, font)
}


#' Internal: order and filter subgroups
#'
#' @noRd
order_subgroups <- function(data, subgroup = FALSE, sort_subgroup_by = "alphabetical") {
  if (!isTRUE(subgroup)) return(data)
  subgroup_order <- setdiff(get_subgroup_order(data, sort_subgroup_by), "Overall")
  data |>
    dplyr::filter(.data$Subgroup %in% subgroup_order) |>
    dplyr::mutate(Subgroup = factor(.data$Subgroup, levels = subgroup_order)) |>
    dplyr::arrange(.data$Subgroup)
}
