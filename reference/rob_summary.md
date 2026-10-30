# Create a Weighted Risk of Bias Summary Plot

Stacked horizontal bar plot showing, for each domain, the proportion of
studies (or of the total weight) at each risk of bias judgement. Uses
the same risk of bias internals as the robbR package.

## Usage

``` r
rob_summary(
  data,
  weight = NULL,
  rob_tool = "rob2",
  show_overall = TRUE,
  derive_overall = FALSE,
  show_na = FALSE,
  wrap = 50,
  font = NULL,
  title = NULL,
  title_align = c("left", "center", "centre", "right"),
  subtitle = NULL
)
```

## Arguments

- data:

  A data frame of domain judgements (see
  [`rob_plot()`](https://blmoran.github.io/bayesfoRest/reference/rob_plot.md)).

- weight:

  \<[`data-masking`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)\>
  Optional weight for each study (e.g. meta-analysis weight or sample
  size). If `NULL` (default), each study counts equally.

- rob_tool:

  Risk of bias tool: `"rob2"` (default), `"robins_i"`, `"robins_e"`,
  `"quadas2"`, `"clarity"`, `"rob_nma"`, `"probast"`, or `"probast_ai"`.

- show_overall:

  Logical. Show the overall bar(s)? (default: `TRUE`).

- derive_overall:

  Logical. Always re-derive overall judgements from the domains?
  (default: `FALSE`).

- show_na:

  Logical. Include "Not applicable" and "Outcome not measured" in the
  bars? If `FALSE` (default) they are excluded from the denominator.

- wrap:

  Width at which domain labels are wrapped.

- font:

  Font family.

- title, subtitle:

  Plot title and subtitle.

- title_align:

  Title alignment: `"left"` (default), `"center"`/`"centre"`, or
  `"right"`.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
rob_summary(rob_data)
rob_summary(rob_data, weight = weight, title = "Risk of bias (weighted)")
rob_summary(robins_data, weight = n, rob_tool = "robins_i")
} # }
```
