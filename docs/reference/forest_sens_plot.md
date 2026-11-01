# Create a Forest Plot with Multi-Prior Pooled-Effect Rows

A classic forest plot – individual study rows exactly like
[`bayes_forest()`](https://blmoran.github.io/bayesfoRest/reference/bayes_forest.md)
– whose "Pooled Effect" section shows one row per prior specification
instead of a single pooled estimate, plus an (optional)
prediction-interval row. Useful for showing how sensitive the pooled
estimate is to the choice of prior directly on the same plot as the
underlying studies. For a fuller sensitivity analysis comparing common-
vs random-effects models, RoB exclusion, or RoBMA bias adjustment (but
without individual study rows), see
[`sensitivity_plot()`](https://blmoran.github.io/bayesfoRest/reference/sensitivity_plot.md).

## Usage

``` r
forest_sens_plot(
  model,
  data,
  measure,
  studyvar = NULL,
  yearvar = NULL,
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
  sort_studies_by = "author",
  priors,
  prior_order = NULL,
  primary_prior = NULL,
  show_prediction = TRUE,
  study_output = c("density", "boxplot", "sq_whisker"),
  pooled_output = c("diamond", "density", "boxplot", "sq_whisker"),
  pred_output = c("diamond", "density", "boxplot", "sq_whisker"),
  label_control = "Control",
  label_intervention = "Intervention",
  title = NULL,
  subtitle = NULL,
  title_align = "left",
  xlim = NULL,
  x_breaks = NULL,
  x_label = NULL,
  null_value = NULL,
  add_rope = FALSE,
  rope_value = NULL,
  rope_color = "grey50",
  color_palette = NULL,
  color_study_posterior = "dodgerblue",
  color_study_posterior_outline = "blue",
  color_pooled_posterior = "blue",
  color_shrinkage_pointinterval = "purple",
  color_shrinkage_outline = "purple",
  color_shrinkage_fill = NULL,
  split_color_by_null = FALSE,
  color_favours_control = "firebrick",
  color_favours_intervention = "dodgerblue",
  add_arm_labels = TRUE,
  reverse_arm_labels = FALSE,
  color_pred_posterior = "orange",
  color_pred_outline = "darkorange2",
  plot_width = 4,
  font = NULL
)
```

## Arguments

- model:

  A fitted Bayesian meta-analysis model object. The study rows support
  the same backends as
  [`bayes_forest()`](https://blmoran.github.io/bayesfoRest/reference/bayes_forest.md):
  `brmsfit`, `meta_stan`, `bayesmeta`, `stanreg`, `stanfit`,
  `CmdStanMCMC`, `bayesma`, `brma`, `rjags`. Refitting under each entry
  in `priors`, however, only supports `brmsfit`, `bayesmeta`,
  `meta_stan` and `brma` (a single
  [`RoBMA::brma()`](https://fbartos.github.io/RoBMA/reference/brma.html)
  fit, not a
  [`RoBMA::RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.html)
  model-averaged ensemble) – other backends raise an error naming the
  unsupported model class. See
  [`bf_priors`](https://blmoran.github.io/bayesfoRest/reference/bf_priors.md)
  for which prior families each of those four backends can fit exactly.

- data:

  A data frame containing the study data used for the meta-analysis.

- measure:

  Character string specifying the effect measure. Must be one of: "OR",
  "HR", "RR", "IRR", "MD", or "SMD".

- studyvar:

  Column name containing study identifiers/authors.

- yearvar:

  Column name containing publication years. Default is NULL.

- n_ctrl:

  Column name containing control group sample sizes.

- n_int:

  Column name containing intervention group sample sizes.

- event_ctrl:

  Column name containing control group event counts.

- event_int:

  Column name containing intervention group event counts.

- mean_ctrl:

  Column name containing control group means.

- mean_int:

  Column name containing intervention group means.

- sd_ctrl:

  Column name containing control group standard deviations.

- sd_int:

  Column name containing intervention group standard deviations.

- time_ctrl:

  Column name containing control group time periods.

- time_int:

  Column name containing intervention group time periods.

- sort_studies_by:

  Character string specifying how to sort studies. Options: "author"
  (default), "year", or "effect".

- priors:

  A named list of prior specifications. Each element should be a list
  with `name` (display label), `mu_prior` (a `bf_prior` object), and
  optionally `tau_prior` (a `bf_prior` object). See
  [`bf_normal`](https://blmoran.github.io/bayesfoRest/reference/bf_priors.md)
  for constructors. Each prior is refit via the same machinery
  [`sensitivity_plot()`](https://blmoran.github.io/bayesfoRest/reference/sensitivity_plot.md)
  uses, and shown as its own pooled-effect row labelled "`name` Prior".

- prior_order:

  Optional character vector of prior keys specifying display order.
  Defaults to the order in `priors`.

- primary_prior:

  Optional character key into `priors`. That prior's row label gets "
  (primary)" appended – use it to flag which prior specification matches
  the prior `model` was actually fit with.

- show_prediction:

  Logical. Show a "Prediction interval" row below the priors (derived
  from `model`'s own pooled estimate and tau). Automatically dropped
  (with a message) for common-effect models, which have no tau. Default
  TRUE.

- study_output:

  Character string specifying how per-study rows are drawn: "sq_whisker"
  (default, matching a classic forest plot), "density", or "boxplot".

- pooled_output:

  Character string specifying how each prior's row is drawn: "diamond"
  (default), "density", "boxplot", or "sq_whisker".

- pred_output:

  Character string specifying how the prediction-interval row is drawn:
  "diamond" (default), "density", "boxplot", or "sq_whisker".

- label_control:

  Character string for control group label. Default "Control".

- label_intervention:

  Character string for intervention group label. Default "Intervention".

- title:

  Character string for the plot title. Default NULL.

- subtitle:

  Character string for the plot subtitle. Default NULL.

- title_align:

  Character string specifying title alignment.

- xlim:

  Numeric vector of length 2 specifying x-axis limits.

- x_breaks:

  Numeric vector specifying custom x-axis break points.

- x_label:

  Character string for the x-axis label. Default NULL (measure-specific
  label).

- null_value:

  Numeric value specifying x-axis value of the null line.

- add_rope:

  Logical indicating whether to add a ROPE. Default FALSE.

- rope_value:

  Numeric vector of length 2 specifying ROPE range.

- rope_color:

  Color for ROPE shading.

- color_palette:

  Character vector of colors for the plot.

- color_study_posterior:

  Color for study posterior densities.

- color_study_posterior_outline:

  Color for study posterior outlines.

- color_pooled_posterior:

  Fill color for each prior's row.

- color_shrinkage_pointinterval:

  Color for shrinkage point intervals.

- color_shrinkage_outline:

  Color for shrinkage plot outlines.

- color_shrinkage_fill:

  Color for shrinkage plot fill.

- split_color_by_null:

  Logical. Split posterior densities by null value.

- color_favours_control:

  Colour for regions favouring control.

- color_favours_intervention:

  Colour for regions favouring intervention.

- add_arm_labels:

  Logical. Show "Favours Control/Intervention" labels.

- reverse_arm_labels:

  Logical. Swap the arm label positions.

- color_pred_posterior:

  Color for the prediction-interval row fill.

- color_pred_outline:

  Color for the prediction-interval row outline.

- plot_width:

  Numeric value specifying the relative width of the plot component.
  Default 4.

- font:

  Character string specifying the font family.

## Value

A patchwork object containing the complete forest plot.

## Examples

``` r
if (FALSE) { # \dontrun{
priors <- list(
  vague = list(
    name      = "Vague",
    mu_prior  = bf_normal(0, 10),
    tau_prior = bf_half_cauchy(0, 1)
  ),
  weakreg = list(
    name      = "Weakly Regularising",
    mu_prior  = bf_normal(0, 1),
    tau_prior = bf_half_cauchy(0, 0.5)
  ),
  informative = list(
    name      = "Informative",
    mu_prior  = bf_normal(0, 0.5),
    tau_prior = bf_half_cauchy(0, 0.3)
  )
)

forest_sens_plot(
  model         = fit,
  data          = dat,
  measure       = "OR",
  studyvar      = Author,
  n_ctrl        = N_Control,
  n_int         = N_Intervention,
  event_ctrl    = Event_Control,
  event_int     = Event_Intervention,
  priors        = priors,
  primary_prior = "weakreg"
)
} # }
```
