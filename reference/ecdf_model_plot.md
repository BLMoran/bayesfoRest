# ECDF Plot Comparing Model Strategies (fixed prior)

ECDF (+ optional density) comparing sensitivity model strategies – All
studies, Excluding High RoB, PET-PEESE, Mixture Model, BMA – for a
single prior. Colour represents the model strategy/section.

## Usage

``` r
ecdf_model_plot(
  model,
  data,
  prior,
  measure,
  study_var = NULL,
  rob_var = NULL,
  exclude_high_rob = FALSE,
  incl_common_effect = FALSE,
  incl_random_effect = TRUE,
  incl_pet_peese = FALSE,
  incl_selection_weight = FALSE,
  incl_robma = FALSE,
  incl_robma_bma = FALSE,
  robma_fit = NULL,
  pet_peese_direction = "negative",
  pet_peese_threshold = 0.1,
  prob_reference = "null",
  null_value = NULL,
  null_range = NULL,
  add_null_range = FALSE,
  color_null_range = "#77bb41",
  label_control = "Control",
  label_intervention = "Intervention",
  title = NULL,
  subtitle = NULL,
  x_label = NULL,
  xlim = NULL,
  x_breaks = NULL,
  color_palette = NULL,
  show_density = TRUE,
  font = NULL
)
```

## Arguments

- model:

  A fitted `brmsfit` meta-analysis model.

- data:

  Study-level data used to fit `model`.

- prior:

  A single prior spec: a list with `name`, `mu_prior`, and optionally
  `tau_prior`. To compare across multiple priors instead of model
  strategies, use
  [`ecdf_prior_plot`](https://blmoran.github.io/bayesfoRest/reference/ecdf_prior_plot.md).

- measure:

  Effect measure ("OR","RR","HR","IRR","MD","SMD").

- study_var, rob_var:

  Optional unquoted column names (study id / RoB).

- exclude_high_rob, incl_common_effect, incl_random_effect,
  incl_pet_peese, incl_selection_weight, incl_robma, incl_robma_bma:

  Logical section toggles. `incl_robma` adds the RoBMA conditional row
  (pooled effect given an effect exists); `incl_robma_bma` adds the
  model-averaged row (can show a point mass at the null) – the two are
  independent.

- robma_fit:

  A single fitted RoBMA object for this prior (required if `incl_robma`,
  `incl_robma_bma`, `incl_pet_peese`, or `incl_selection_weight`).

- pet_peese_direction, pet_peese_threshold:

  PET-PEESE options.

- prob_reference:

  "null" or "null_range".

- null_value, null_range, add_null_range, color_null_range:

  Null / ROPE options.

- label_control, label_intervention:

  Arm labels for the density panel.

- title, subtitle, x_label, xlim, x_breaks, color_palette, show_density,
  font:

  Display options.

## Value

A ggplot (or patchwork) object.
