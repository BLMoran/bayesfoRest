# Funnel Plot (multi-backend)

Publication-ready funnel plot for meta-analysis, displaying study-level
effect sizes against their precision (standard error), styled to match
the bayesfoRest package. The plot itself only needs study-level
`yi`/`sei` and a single pooled estimate, so it accepts any fitted model
backend that bayesfoRest supports for the pooled line: `brmsfit`,
`bayesmeta`, MetaStan / `stanfit` (rstan), `CmdStanMCMC` (cmdstanr),
`rjags`/`runjags`, or `bayesma`. The pooled value is obtained via the
same backend-agnostic accessor used by
[`bayes_forest()`](https://blmoran.github.io/bayesfoRest/reference/bayes_forest.md).

## Usage

``` r
funnel_plot(
  model,
  data,
  measure,
  studyvar = NULL,
  yearvar = NULL,
  yi_col = NULL,
  sei_col = NULL,
  n_ctrl = NULL,
  n_int = NULL,
  event_ctrl = NULL,
  event_int = NULL,
  mean_ctrl = NULL,
  mean_int = NULL,
  sd_ctrl = NULL,
  sd_int = NULL,
  time_ctrl = NULL,
  time_int = NULL,
  subgroup = FALSE,
  subgroupvar = NULL,
  contour_lines = TRUE,
  contour_alpha = 0.15,
  pooled_line = c("posterior", "fixed", "none"),
  color_points = "dodgerblue",
  color_points_outline = "blue",
  color_pooled = "blue",
  color_null = "black",
  color_contour = "grey50",
  point_size = 3,
  xlim = NULL,
  x_breaks = NULL,
  null_value = NULL,
  add_null_line = FALSE,
  label_studies = FALSE,
  label_size = 3,
  title = NULL,
  subtitle = NULL,
  x_label = NULL,
  title_align = "left",
  font = NULL
)
```

## Arguments

- model:

  A fitted meta-analysis model (any supported backend). Used only to
  draw the pooled-estimate line (when `pooled_line = "posterior"`).

- data:

  A data frame of study-level data. Must either contain the columns
  named by `yi_col`/`sei_col`, or the arm columns needed to compute
  `yi`/`sei` for `measure` (as in
  [`bayes_forest()`](https://blmoran.github.io/bayesfoRest/reference/bayes_forest.md)).

- measure:

  Character; effect measure. One of `"OR"`, `"RR"`, `"HR"`, `"IRR"`,
  `"MD"`, `"SMD"`.

- studyvar:

  Unquoted column with study identifiers.

- yearvar:

  Unquoted column with study year (optional, used in labels).

- yi_col, sei_col:

  Optional character names of pre-computed effect and SE columns (e.g.
  difference-of-medians `yi`/`sei`). If supplied, `yi`/`sei` are taken
  directly and no arm columns are needed.

- n_ctrl, n_int, event_ctrl, event_int, mean_ctrl, mean_int, sd_ctrl,
  sd_int, time_ctrl, time_int:

  Unquoted arm columns, used only when `yi_col`/`sei_col` are not
  supplied.

- subgroup:

  Logical; colour points by subgroup.

- subgroupvar:

  Unquoted subgroup column (when `subgroup = TRUE`).

- contour_lines:

  Logical; add 90/95/99% significance contour shading.

- contour_alpha:

  Numeric; base alpha for contour shading.

- pooled_line:

  One of `"posterior"` (posterior median of the pooled effect,
  backend-agnostic), `"fixed"` (inverse-variance fixed-effect estimate
  from the data), or `"none"`.

- color_points, color_points_outline, color_pooled, color_null,
  color_contour:

  Colours.

- point_size, label_size:

  Sizes.

- xlim, x_breaks, x_label:

  Axis limits / breaks / label.

- null_value, add_null_line:

  Null-effect line.

- label_studies:

  Logical; label points with study (year).

- title, subtitle, title_align, font:

  Display options.

## Value

A ggplot object.
