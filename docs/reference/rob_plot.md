# Create Risk of Bias Plot

Creates a publication-ready risk of bias plot using gt tables, for
single studies or subgroup analyses. Uses the same risk of bias
internals as the robbR package.

## Usage

``` r
rob_plot(
  data,
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
  subtitle = NULL
)
```

## Arguments

- data:

  A data frame with a study, year and domain columns for the chosen
  `rob_tool` (see Details). Judgements are standardised (e.g. "low risk"
  -\> "Low"); `NA` -\> "Not applicable". Judgements not recognised by
  the tool are shown as blank cells, with a warning.

- studyvar, yearvar:

  \<[`data-masking`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)\>
  Study label and year columns.

- subgroupvar:

  \<[`data-masking`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)\>
  Subgroup column, used when `subgroup = TRUE`.

- effect:

  \<[`data-masking`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)\>
  Effect size column, used when sorting by `"effect"`.

- rob_tool:

  Risk of bias tool: `"rob2"` (default), `"robins_i"`, `"robins_e"`,
  `"quadas2"`, `"clarity"`, `"rob_nma"`, `"probast"`, or `"probast_ai"`.

- sort_studies_by:

  How to sort studies: `"author"` (default), `"year"`, or `"effect"`.

- subgroup:

  Logical. Create a subgroup plot? (default: `FALSE`).

- sort_subgroup_by:

  Subgroup ordering: `"alphabetical"` (default), `"effect"`, or a
  character vector of subgroup names in the desired order.

- show_overall:

  Logical. Show the overall column(s)? (default: `TRUE`). Overall
  columns absent from `data` are derived from the domains (worst
  judgement).

- derive_overall:

  Logical. Always re-derive overall column(s) from the domains,
  overwriting supplied values? (default: `FALSE`).

- add_rob_legend:

  Logical. Add a legend? (default: `FALSE`).

- study_label:

  Header for the study column. Defaults to `"Trial"` for randomised
  trial tools (`"rob2"`, `"clarity"`), `"Review"` for `"rob_nma"`, and
  `"Study"` otherwise.

- font:

  Font family.

- title, subtitle:

  Plot title and subtitle.

- title_align:

  Title alignment: `"left"` (default), `"center"`/`"centre"`, or
  `"right"`.

## Value

A patchwork composition of the risk of bias table (and legend).

## Details

Domain columns by tool (optional columns in brackets are shown only if
present; `NA` rows are shown as "Not applicable"):

- rob2:

  D1, (D1b), D2-D5, (DS), Overall.

- robins_i:

  ROBINS-I V2: D1-D6, Overall.

- robins_e:

  D1-D7, Overall.

- quadas2:

  D1-D4, Overall.

- clarity:

  D1-D10, (D11), (D12). No overall judgement.

- rob_nma:

  (R1-R3), D1-D3, Results, Conclusions.

- probast, probast_ai:

  D1-D4, (A1-A3), Overall, Applicability.

## Examples

``` r
if (FALSE) { # \dontrun{
rob_data <- data.frame(
  Author = c("Smith", "Jones", "Brown"),
  Year = c(2020, 2021, 2022),
  D1 = c("Low", "High", "Some concerns"),
  D2 = c("Low", "Low", "High"),
  D3 = c("Some concerns", "Low", "Low"),
  D4 = "Low",
  D5 = "Low")

rob_plot(rob_data, add_rob_legend = TRUE)
} # }
```
