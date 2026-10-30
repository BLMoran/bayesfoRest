# Risk of bias internals, vendored from the robbR package
# (https://github.com/BLMoran/robbR). Keep in sync with robbR/R/rob_spec.R,
# rob_tools.R, rob_judgement.R and rob_table_fn.R.

#' Risk of bias tool specification
#'
#' S7 class describing a risk of bias tool: its domains, judgement levels, how
#' judgements are standardised, and how overall judgements are derived.
#'
#' @noRd
rob_spec <- S7::new_class(
  "rob_spec",
  properties = list(
    name = S7::class_character,
    label = S7::class_character,
    domains = S7::class_character,
    levels = S7::class_character,
    severity = S7::class_integer,
    synonyms = S7::class_data.frame,
    optional = S7::class_character,
    always = S7::class_character,
    overall = S7::class_character,
    overall_from = S7::class_list,
    overall_map = S7::class_character,
    overall_synonyms = S7::class_data.frame,
    spanners = S7::class_list,
    study_label = S7::class_character,
    required = S7::new_property(
      S7::class_character,
      getter = \(self) setdiff(names(self@domains), self@optional))),
  constructor = function(name, label, domains, levels, severity, synonyms, optional = character(), always = character(), overall = c(Overall = "Overall risk of bias"), overall_from = NULL, overall_map = NULL, overall_synonyms = NULL, spanners = list(), study_label = "Study") {
    levels <- c(levels, "Not applicable" = "Not applicable", ONM = "Outcome not measured")
    S7::new_object(
      S7::S7_object(),
      name = name,
      label = label,
      domains = domains,
      levels = levels,
      severity = rlang::set_names(as.integer(severity), names(severity)),
      synonyms = rob_synonyms(synonyms, names(levels)),
      optional = optional,
      always = always,
      overall = overall,
      overall_from = overall_from %||% list(Overall = names(domains)),
      overall_map = overall_map %||% rlang::set_names(names(levels)),
      overall_synonyms = rob_synonyms(overall_synonyms %||% synonyms, names(levels)),
      spanners = spanners,
      study_label = study_label)
  },
  validator = function(self) {
    c(
      if (!all(self@optional %in% names(self@domains))) "@optional must be domain names.",
      if (!all(self@always %in% self@optional)) "@always must be optional domains.",
      if (!all(names(self@overall_from) %in% names(self@overall))) "@overall_from must be named by overall columns.",
      if (!all(names(self@severity) %in% names(self@levels))) "@severity must be named by levels in @levels.")})


#' Coerce to a risk of bias tool specification
#'
#' @noRd
as_rob_spec <- S7::new_generic("as_rob_spec", "x")

S7::method(as_rob_spec, rob_spec) <- function(x) x

S7::method(as_rob_spec, S7::class_character) <- function(x) get_rob_spec(x)


#' Overall judgement for one set of domain judgements
#'
#' @noRd
rob_overall_judgement <- S7::new_generic("rob_overall_judgement", "spec")

S7::method(rob_overall_judgement, rob_spec) <- function(spec, x) {
  severity <- tibble::tibble(level = names(spec@severity), severity = unname(spec@severity), rank = seq_along(spec@severity))
  judged <- tibble::tibble(level = dplyr::coalesce(unname(spec@overall_map[x]), x)) |>
    dplyr::inner_join(severity, by = "level", relationship = "many-to-one") |>
    dplyr::slice_max(.data$severity, n = 1, with_ties = TRUE) |>
    dplyr::slice_min(.data$rank, n = 1, with_ties = FALSE)
  if (nrow(judged) == 0) "Not applicable" else judged$level
}


#' Internal: synonym lookup table
#'
#' @noRd
rob_synonyms <- function(synonyms, levels) {
  c(synonyms, list("Not applicable" = c("na", "n/a", "not applicable", ""), ONM = c("onm", "outcome not measured"))) |>
    purrr::imap(\(synonym, level) tibble::tibble(synonym = synonym, level = level)) |>
    purrr::list_rbind() |>
    dplyr::bind_rows(tibble::tibble(synonym = stringr::str_to_lower(levels), level = levels)) |>
    dplyr::slice_head(n = 1, by = "synonym")
}


#' Built-in risk of bias tools
#'
#' @noRd
rob_tools <- function() {
  c("rob2", "robins_i", "robins_e", "quadas2", "clarity", "rob_nma", "probast", "probast_ai")
}


#' Get a built-in risk of bias tool specification
#'
#' @noRd
get_rob_spec <- function(rob_tool = "rob2") {
  rob_tool <- rlang::arg_match(rob_tool, rob_tools())
  switch(
    rob_tool,
    rob2 = spec_rob2(),
    robins_i = spec_robins_i(),
    robins_e = spec_robins_e(),
    quadas2 = spec_quadas2(),
    clarity = spec_clarity(),
    rob_nma = spec_rob_nma(),
    probast = spec_probast("probast", "PROBAST", "Participants", "Overall risk of bias"),
    probast_ai = spec_probast("probast_ai", "PROBAST+AI", "Participants and data sources", "Overall quality (development) / risk of bias (evaluation)"))
}


syn_low <- c("low", "low risk", "low risk of bias", "low concern", "low concerns")
syn_some <- c("some concerns", "some concern", "concerns")
syn_high <- c("high", "high risk", "high risk of bias", "high concern", "high concerns")
syn_unclear <- c("unclear", "unclear risk", "unclear risk of bias", "no information", "ni")
syn_low_except <- c(
  "low (except for concerns about uncontrolled confounding)",
  "low risk of bias (except for concerns about uncontrolled confounding)",
  "low risk of bias except for concerns about uncontrolled confounding",
  "low except for concerns about uncontrolled confounding",
  "low (except for uncontrolled confounding)",
  "low (except confounding)",
  "low except confounding")

robins_low_levels <- c(
  "Low" = "Low risk of bias",
  "Low (except confounding)" = "Low risk of bias (except for concerns about uncontrolled confounding)")


spec_rob2 <- function() {
  rob_spec(
    name = "rob2",
    study_label = "Trial",
    label = "RoB 2",
    domains = c(
      D1 = "Bias arising from the randomization process",
      D1b = "Bias arising from the timing of identification or recruitment of participants (cluster-randomized trials)",
      D2 = "Bias due to deviations from intended interventions",
      D3 = "Bias due to missing outcome data",
      D4 = "Bias in measurement of the outcome",
      D5 = "Bias in selection of the reported result",
      DS = "Bias arising from period and carryover effects (crossover trials)"),
    optional = c("D1b", "DS"),
    levels = c(Low = "Low risk of bias", "Some concerns" = "Some concerns", High = "High risk of bias"),
    severity = c(Low = 1L, "Some concerns" = 2L, High = 3L),
    synonyms = list(Low = syn_low, "Some concerns" = syn_some, High = syn_high))
}


spec_robins_i <- function() {
  rob_spec(
    name = "robins_i",
    label = "ROBINS-I V2",
    domains = c(
      D1 = "Bias due to confounding",
      D2 = "Bias in classification of interventions",
      D3 = "Bias in selection of participants into the study (or into the analysis)",
      D4 = "Bias due to missing data",
      D5 = "Bias in measurement of the outcome",
      D6 = "Bias in selection of the reported result"),
    levels = c(
      robins_low_levels,
      Moderate = "Moderate risk of bias",
      Serious = "Serious risk of bias",
      Critical = "Critical risk of bias"),
    severity = c(Low = 1L, "Low (except confounding)" = 2L, Moderate = 3L, Serious = 4L, Critical = 5L),
    synonyms = list(
      Low = syn_low,
      "Low (except confounding)" = syn_low_except,
      Moderate = c("moderate", "moderate risk", "moderate risk of bias"),
      Serious = c("serious", "serious risk", "serious risk of bias"),
      Critical = c("critical", "critical risk", "critical risk of bias")))
}


spec_robins_e <- function() {
  rob_spec(
    name = "robins_e",
    label = "ROBINS-E",
    domains = c(
      D1 = "Bias due to confounding",
      D2 = "Bias arising from measurement of the exposure",
      D3 = "Bias in selection of participants into the study (or into the analysis)",
      D4 = "Bias due to post-exposure interventions",
      D5 = "Bias due to missing data",
      D6 = "Bias arising from measurement of the outcome",
      D7 = "Bias in selection of the reported result"),
    levels = c(
      robins_low_levels,
      "Some concerns" = "Some concerns",
      High = "High risk of bias",
      "Very high" = "Very high risk of bias"),
    severity = c(Low = 1L, "Low (except confounding)" = 2L, "Some concerns" = 3L, High = 4L, "Very high" = 5L),
    synonyms = list(
      Low = syn_low,
      "Low (except confounding)" = syn_low_except,
      "Some concerns" = syn_some,
      High = syn_high,
      "Very high" = c("very high", "very high risk", "very high risk of bias")))
}


spec_quadas2 <- function() {
  rob_spec(
    name = "quadas2",
    label = "QUADAS-2",
    domains = c(
      D1 = "Patient selection",
      D2 = "Index test(s)",
      D3 = "Reference standard",
      D4 = "Flow and timing"),
    levels = c(Low = "Low risk of bias", Unclear = "Unclear risk of bias", High = "High risk of bias"),
    severity = c(Low = 1L, Unclear = 2L, High = 3L),
    synonyms = list(Low = syn_low, Unclear = syn_unclear, High = syn_high))
}


spec_clarity <- function() {
  rob_spec(
    name = "clarity",
    study_label = "Trial",
    label = "CLARITY",
    domains = c(
      D1 = "Allocation sequence adequately generated",
      D2 = "Allocation adequately concealed",
      D3 = "Patients blinded",
      D4 = "Healthcare providers blinded",
      D5 = "Data collectors blinded",
      D6 = "Outcome assessors blinded",
      D7 = "Data analysts blinded",
      D8 = "Missing outcome data infrequent",
      D9 = "No selective outcome reporting",
      D10 = "Free of other bias",
      D11 = "No selective enrolment within clusters (cluster trials)",
      D12 = "No changes between intervention and control periods (cluster cross-over trials)"),
    optional = c("D11", "D12"),
    always = c("D11", "D12"),
    overall = character(),
    overall_from = list(),
    levels = c(
      "Definitely low" = "Definitely low risk",
      "Probably low" = "Probably low risk",
      "Probably high" = "Probably high risk",
      "Definitely high" = "Definitely high risk"),
    severity = c("Definitely low" = 1L, "Probably low" = 2L, "Probably high" = 3L, "Definitely high" = 4L),
    synonyms = list(
      "Definitely low" = c("low", "definitely low"),
      "Probably low" = "probably low",
      "Probably high" = "probably high",
      "Definitely high" = c("high", "definitely high")))
}


spec_rob_nma <- function() {
  rob_spec(
    name = "rob_nma",
    study_label = "Review",
    label = "RoB NMA",
    domains = c(
      R1 = "ROBIS: Study eligibility criteria",
      R2 = "ROBIS: Identification and selection of studies",
      R3 = "ROBIS: Data collection and study appraisal",
      D1 = "Interventions and network geometry",
      D2 = "Effect modifiers",
      D3 = "Statistical synthesis"),
    optional = c("R1", "R2", "R3"),
    overall = c(Results = "Overall risk of bias in the NMA results", Conclusions = "Concerns about the authors' conclusions"),
    overall_from = list(Results = c("R1", "R2", "R3", "D1", "D2", "D3"), Conclusions = NULL),
    levels = c(Low = "Low risk of bias / no concerns", "Some concerns" = "Some concerns", Unclear = "Unclear", High = "High risk of bias"),
    severity = c(Low = 1L, "Some concerns" = 2L, Unclear = 2L, High = 3L),
    synonyms = list(Low = c(syn_low, "no concerns", "no concern"), "Some concerns" = syn_some, Unclear = syn_unclear, High = syn_high),
    spanners = list(ROBIS = c("R1", "R2", "R3"), "RoB NMA" = c("D1", "D2", "D3"), Overall = c("Results", "Conclusions")))
}


spec_probast <- function(name, label, participants, overall_label) {
  spanner_label <- if (name == "probast_ai") "Quality / risk of bias" else "Risk of bias"
  rob_spec(
    name = name,
    label = label,
    domains = c(
      D1 = participants,
      D2 = "Predictors",
      D3 = "Outcome",
      D4 = "Analysis",
      A1 = stringr::str_c("Applicability: ", stringr::str_to_lower(participants)),
      A2 = "Applicability: predictors",
      A3 = "Applicability: outcome"),
    optional = c("A1", "A2", "A3"),
    overall = c(Overall = overall_label, Applicability = "Overall applicability"),
    overall_from = list(Overall = c("D1", "D2", "D3", "D4"), Applicability = c("A1", "A2", "A3")),
    levels = c(Low = "Low risk / concern", Unclear = "Unclear", High = "High risk / concern"),
    severity = c(Low = 1L, Unclear = 2L, High = 3L),
    synonyms = list(Low = syn_low, Unclear = syn_unclear, High = syn_high),
    spanners = rlang::list2(
      !!spanner_label := c("D1", "D2", "D3", "D4", "Overall"),
      Applicability = c("A1", "A2", "A3", "Applicability")))
}

#' Risk of bias columns present in a data frame
#'
#' @noRd
rob_columns <- function(data, rob_tool = "rob2") {
  spec <- as_rob_spec(rob_tool)
  data |>
    dplyr::select(dplyr::any_of(c(names(spec@domains), names(spec@overall)))) |>
    names()
}


#' Internal: standardise judgements against a synonym lookup
#'
#' @noRd
standardise_judgement <- function(x, lookup) {
  tibble::tibble(value = as.character(x)) |>
    dplyr::mutate(synonym = stringr::str_to_lower(stringr::str_squish(.data$value))) |>
    dplyr::left_join(lookup, by = "synonym", relationship = "many-to-one") |>
    dplyr::mutate(level = dplyr::case_when(is.na(.data$value) ~ "Not applicable", .default = dplyr::coalesce(.data$level, .data$value))) |>
    dplyr::pull("level")
}


#' Internal: validate, standardise and derive overall columns
#'
#' @noRd
prepare_rob_data <- function(data, spec, derive_overall = FALSE) {
  missing_cols <- setdiff(spec@required, names(data))
  if (length(missing_cols) > 0) cli::cli_abort("Missing {spec@label} column{?s}: {.field {missing_cols}}.")
  data <- data |>
    dplyr::mutate(!!!rlang::rep_named(setdiff(spec@always, names(data)), list(NA_character_))) |>
    dplyr::mutate(dplyr::across(dplyr::any_of(names(spec@domains)), \(x) standardise_judgement(x, spec@synonyms))) |>
    purrr::reduce(names(spec@overall), \(df, name) add_overall(df, spec, name, derive_overall), .init = _)
  warn_unknown_judgements(data, spec)
  data
}


#' Internal: derive an overall judgement column from its component domains
#'
#' Applies rob_overall_judgement() row-wise across the given columns.
#'
#' @noRd
rob_overall_column <- function(data, spec, cols, name) {
  data |>
    dplyr::rowwise() |>
    dplyr::mutate("{name}" := rob_overall_judgement(spec, dplyr::c_across(dplyr::all_of(cols)))) |>
    dplyr::ungroup()
}


#' Internal: derive an overall column if absent (or forced), else standardise it
#'
#' @noRd
add_overall <- function(data, spec, name, derive_overall = FALSE) {
  cols <- dplyr::intersect(spec@overall_from[[name]], names(data))
  if ((isTRUE(derive_overall) || !name %in% names(data)) && length(cols) > 0) return(rob_overall_column(data, spec, cols = cols, name = name))
  if (!name %in% names(data)) return(data)
  dplyr::mutate(data, "{name}" := standardise_judgement(.data[[name]], spec@overall_synonyms))
}


#' Internal: warn about judgements that will be shown as blank cells
#'
#' @noRd
warn_unknown_judgements <- function(data, spec) {
  unknown <- data |>
    dplyr::select(dplyr::all_of(rob_columns(data, spec))) |>
    purrr::map(as.character) |>
    purrr::list_c() |>
    unique() |>
    setdiff(names(spec@levels))
  if (length(unknown) > 0) cli::cli_warn("Unrecognised {spec@label} judgement{?s} shown as blank cell{?s}: {.val {unknown}}.")
  invisible(unknown)
}


#' Internal: gt table of risk of bias judgements
#'
#' @noRd
rob_table_fn <- function(df, spec, rob_cols, subgroup = FALSE, study_label = NULL, font = NULL) {
  study_label <- study_label %||% spec@study_label
  df <- dplyr::mutate(df, Author = dplyr::if_else(is.na(.data$Year), as.character(.data$Author), stringr::str_c(.data$Author, " (", .data$Year, ")")))
  base_table <- if (isTRUE(subgroup)) {
    gt::gt(df, groupname_col = "Subgroup") |>
      gt::tab_options(row_group.font.weight = "bold") |>
      gt::tab_style(style = gt::cell_fill(color = "grey95"), locations = gt::cells_row_groups()) |>
      gt::cols_label(Author = gt::md(stringr::str_c("Subgroup/ \n", study_label))) |>
      gt::tab_style(style = gt::cell_text(indent = gt::px(10)), locations = gt::cells_body(columns = "Author"))
  } else {
    gt::gt(df) |>
      gt::cols_label(Author = study_label)
  }
  base_table |>
    gt::cols_align(align = "left", columns = "Author") |>
    gt::cols_align(align = "center", columns = dplyr::all_of(rob_cols)) |>
    gt::tab_options(column_labels.font.weight = "bold", table.border.top.color = "white", table.border.bottom.color = "white", table.font.names = font) |>
    gt::cols_hide(columns = "Year") |>
    gt::opt_table_lines(extent = "none") |>
    gt::text_transform(locations = gt::cells_body(columns = dplyr::all_of(rob_cols)), fn = \(x) rob_icon(x, levels = names(spec@levels), height = "1.1em")) |>
    add_rob_spanners(spec, rob_cols)
}


#' Internal: column spanners (shown when at least two groups are present)
#'
#' @noRd
add_rob_spanners <- function(table, spec, rob_cols) {
  spanners <- spec@spanners |>
    purrr::map(\(x) dplyr::intersect(x, rob_cols)) |>
    purrr::keep(\(x) length(x) > 0)
  if (length(spanners) < 2) return(table)
  purrr::reduce2(names(spanners), spanners, \(tbl, label, cols) gt::tab_spanner(tbl, label = label, columns = dplyr::all_of(cols), id = stringr::str_c("spanner_", label)), .init = table) |>
    gt::tab_style(style = gt::cell_text(weight = "bold"), locations = gt::cells_column_spanners())
}


#' Internal: risk of bias legend
#'
#' Only the domains and judgement levels present in the plot are listed.
#'
#' @noRd
rob_legend_fn <- function(spec, domains_present = NULL, levels_present = NULL, font = NULL) {
  domains <- spec@domains[names(spec@domains) %in% (domains_present %||% names(spec@domains))]
  levels <- spec@levels[names(spec@levels) %in% (levels_present %||% names(spec@severity))]
  legend <- dplyr::bind_rows(
    tibble::tibble(Code = stringr::str_c("(", names(domains), ")"), Description = unname(domains), type = "domain"),
    tibble::tibble(Code = "", Description = "", type = "spacer"),
    tibble::tibble(Code = names(levels), Description = unname(levels), type = "level")) |>
    gt::gt() |>
    gt::cols_hide(columns = "type") |>
    gt::cols_label(Code = "", Description = "") |>
    gt::tab_header(title = "Risk of bias legend") |>
    gt::tab_style(style = gt::cell_text(weight = "bold", size = gt::px(16)), locations = gt::cells_title(groups = "title")) |>
    gt::tab_style(style = gt::cell_text(weight = "bold", size = gt::px(12)), locations = gt::cells_body(columns = "Code", rows = .data$type == "domain")) |>
    gt::tab_style(style = "padding-left:2px;", locations = gt::cells_body(columns = "Code")) |>
    gt::cols_width(Code ~ gt::px(40), Description ~ gt::px(300)) |>
    gt::tab_options(table.border.top.color = "white", table.border.bottom.color = "white", heading.align = "left", heading.padding = gt::px(6), table.font.size = gt::px(12), table.font.names = font) |>
    gt::opt_table_lines(extent = "none")
  if (length(levels) == 0) return(legend)
  legend |>
    gt::text_transform(locations = gt::cells_body(columns = "Code", rows = .data$type == "level"), fn = \(x) rob_icon(x, levels = names(spec@levels), height = "1.5em")) |>
    gt::tab_style(style = gt::cell_text(weight = "bold", style = "italic"), locations = gt::cells_body(columns = "Description", rows = .data$type == "level"))
}


#' Internal: icon for each judgement level
#'
#' @noRd
rob_icon_table <- function() {
  tibble::tribble(
    ~level, ~icon, ~fill, ~prefer,
    "Low", "circle-minus", "#77bb41", "solid",
    "Low (except confounding)", "circle-minus", "#b5d98f", "solid",
    "Some concerns", "circle-question", "#f5ec00", "solid",
    "Moderate", "circle-question", "#f5ec00", "solid",
    "Unclear", "circle-question", "#5b9bd5", "solid",
    "High", "circle-plus", "#e32400", "solid",
    "Serious", "circle-plus", "#e32400", "solid",
    "Very high", "circle-exclamation", "#8b0000", "solid",
    "Critical", "circle-exclamation", "#8b0000", "solid",
    "Definitely low", "circle-minus", "#3d7a1f", "solid",
    "Probably low", "circle-minus", "#a9d18e", "solid",
    "Probably high", "circle-plus", "#f5a000", "solid",
    "Definitely high", "circle-plus", "#e32400", "solid",
    "Not applicable", "circle", "#d9d9d9", "solid",
    "ONM", "circle", "#9e9e9e", "regular")
}


#' Internal: map judgements to icons
#'
#' Judgements not in `levels` (or missing) are rendered as empty cells.
#'
#' @noRd
rob_icon <- function(x, levels = NULL, height = "1.1em") {
  icons <- rob_icon_table()
  if (!is.null(levels)) icons <- dplyr::filter(icons, .data$level %in% levels)
  tibble::tibble(level = as.character(x)) |>
    dplyr::left_join(icons, by = "level", relationship = "many-to-one") |>
    dplyr::select(-"level") |>
    purrr::pmap_chr(\(icon, fill, prefer) if (is.na(icon)) "" else as.character(fontawesome::fa(icon, fill = fill, height = height, prefer_type = prefer)))
}


#' Internal: add title and subtitle to a patchwork plot
#'
#' @noRd
add_rob_titles <- function(plot, title = NULL, subtitle = NULL, title_align = "left", font = NULL) {
  if (is.null(title) && is.null(subtitle)) return(plot)
  hjust <- dplyr::case_when(title_align == "left" ~ 0, title_align %in% c("center", "centre") ~ 0.5, .default = 1)
  title_theme <- ggplot2::theme(
    plot.title = if (!is.null(title)) ggplot2::element_text(size = 16, face = "bold", hjust = hjust, margin = ggplot2::margin(b = if (is.null(subtitle)) 10 else 5), family = font),
    plot.subtitle = if (!is.null(subtitle)) ggplot2::element_text(size = 14, hjust = hjust, margin = ggplot2::margin(b = 10), color = "gray30", family = font))
  plot + patchwork::plot_annotation(title = title, subtitle = subtitle, theme = title_theme)
}


# ---- bayesfoRest-specific helpers (forest plot RoB columns) ----

#' Judgements treated as "high" risk of bias when `exclude_high_rob = TRUE`
#'
#' @noRd
rob_high_levels <- function() {
  c("High", "Serious", "Critical", "Very high", "Definitely high")
}


#' Internal function to get columns that carry the Risk/of/Bias header text
#'
#' @noRd
get_rob_text_columns <- function(rob_cols) {
  n_rob_cols <- length(rob_cols)
  text_positions <- if (n_rob_cols >= 3) {
    start_pos <- ceiling((n_rob_cols - 3) / 2) + 1
    start_pos:(start_pos + 2)
  } else {
    seq_len(n_rob_cols)
  }
  rob_cols[text_positions]
}


#' Internal function to create column labels for the right-hand forest table
#'
#' @noRd
create_rob_labels <- function(rob_cols, measure, effect_label = NULL, effect_type = "random_effect") {
  effect_label <- if (is.null(effect_label)) measure else effect_label

  base_labels <- list(
    weighted_effect = gt::md(paste("Shrinkage", effect_label, "\n", "[95% CrI]")),
    unweighted_effect = gt::md(paste("Observed", effect_label, "\n", "[95% CrI]"))
  )
  if (effect_type == "common_effect") base_labels$weighted_effect <- NULL

  text_cols <- get_rob_text_columns(rob_cols)
  rob_text <- c("Risk", "of", "Bias")[seq_along(text_cols)]

  rob_labels <- purrr::map(rob_cols, ~ {
    if (.x %in% text_cols) {
      gt::md(paste0("**", rob_text[match(.x, text_cols)], "**\n", .x))
    } else {
      gt::md(paste0("&nbsp;\n", .x))
    }
  }) |>
    purrr::set_names(rob_cols)

  c(base_labels, rob_labels)
}
