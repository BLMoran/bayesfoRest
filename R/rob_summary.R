#' Create a Weighted Risk of Bias Summary Plot
#'
#' Stacked horizontal bar plot showing, for each domain, the proportion of
#' studies (or of the total weight) at each risk of bias judgement. Uses the
#' same risk of bias internals as the robbR package.
#'
#' @param data A data frame of domain judgements (see [rob_plot()]).
#' @param weight <[`data-masking`][dplyr::dplyr_data_masking]> Optional weight
#'   for each study (e.g. meta-analysis weight or sample size). If `NULL`
#'   (default), each study counts equally.
#' @param rob_tool Risk of bias tool: `"rob2"` (default), `"robins_i"`,
#'   `"robins_e"`, `"quadas2"`, `"clarity"`, `"rob_nma"`, `"probast"`, or
#'   `"probast_ai"`.
#' @param show_overall Logical. Show the overall bar(s)? (default: `TRUE`).
#' @param derive_overall Logical. Always re-derive overall judgements from the
#'   domains? (default: `FALSE`).
#' @param show_na Logical. Include "Not applicable" and "Outcome not measured"
#'   in the bars? If `FALSE` (default) they are excluded from the denominator.
#' @param wrap Width at which domain labels are wrapped.
#' @param font Font family.
#' @param title,subtitle Plot title and subtitle.
#' @param title_align Title alignment: `"left"` (default), `"center"`/`"centre"`,
#'   or `"right"`.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' rob_summary(rob_data)
#' rob_summary(rob_data, weight = weight, title = "Risk of bias (weighted)")
#' rob_summary(robins_data, weight = n, rob_tool = "robins_i")
#' }
#'
#' @importFrom rlang %||%
#' @export
rob_summary <- function(data,
                        weight = NULL,
                        rob_tool = "rob2",
                        show_overall = TRUE,
                        derive_overall = FALSE,
                        show_na = FALSE,
                        wrap = 50,
                        font = NULL,
                        title = NULL,
                        title_align = c("left", "center", "centre", "right"),
                        subtitle = NULL) {
  title_align <- rlang::arg_match(title_align)
  spec <- as_rob_spec(rob_tool)
  data <- data |>
    dplyr::mutate(.weight = {{ weight }} %||% 1) |>
    prepare_rob_data(spec, derive_overall = derive_overall)
  if (any(is.na(data$.weight) | data$.weight < 0)) cli::cli_abort("{.arg weight} must be non-missing and non-negative.")
  rob_cols <- rob_columns(data, spec)
  if (isFALSE(show_overall)) rob_cols <- setdiff(rob_cols, names(spec@overall))
  proportions <- rob_proportions(data, spec, rob_cols, show_na = show_na)
  levels_present <- dplyr::intersect(names(spec@levels), as.character(proportions$level))
  fills <- rob_icon_table() |>
    dplyr::select("level", "fill") |>
    tibble::deframe()
  labels <- c(spec@domains, spec@overall)[rob_cols] |>
    stringr::str_wrap(width = wrap) |>
    rlang::set_names(rob_cols)
  plot <- ggplot2::ggplot(proportions, ggplot2::aes(x = .data$proportion, y = .data$domain, fill = .data$level)) +
    ggplot2::geom_col(position = ggplot2::position_fill(reverse = TRUE), colour = "black", linewidth = 0.4, width = 0.7) +
    ggplot2::scale_fill_manual(values = fills, breaks = levels_present, labels = stringr::str_wrap(spec@levels[levels_present], width = 22)) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(), breaks = seq(0, 1, 0.25), expand = ggplot2::expansion(mult = c(0.01, 0.03))) +
    ggplot2::scale_y_discrete(labels = labels) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = ceiling(length(levels_present) / 3), byrow = TRUE)) +
    ggplot2::theme_classic(base_size = 12, base_family = font %||% "") +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "black", size = 12),
      legend.position = "bottom",
      legend.location = "plot",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11),
      legend.key.size = ggplot2::unit(0.8, "cm"),
      legend.key.spacing.x = ggplot2::unit(0.5, "cm"),
      legend.background = ggplot2::element_rect(colour = "black", linewidth = 0.4),
      legend.margin = ggplot2::margin(8, 12, 8, 12))
  add_rob_titles(plot, title, subtitle, title_align, font)
}


#' Internal: weighted proportion of each judgement within each domain
#'
#' @noRd
rob_proportions <- function(data, spec, rob_cols, show_na = FALSE) {
  excluded <- if (isTRUE(show_na)) character() else c("Not applicable", "ONM")
  data |>
    dplyr::select(dplyr::all_of(c(rob_cols, ".weight"))) |>
    tidyr::pivot_longer(dplyr::all_of(rob_cols), names_to = "domain", values_to = "level") |>
    dplyr::filter(.data$level %in% setdiff(names(spec@levels), excluded)) |>
    dplyr::summarise(weight = sum(.data$.weight), .by = c("domain", "level")) |>
    dplyr::mutate(proportion = .data$weight / sum(.data$weight), .by = "domain") |>
    dplyr::mutate(domain = factor(.data$domain, levels = rev(rob_cols)), level = factor(.data$level, levels = names(spec@levels)))
}
