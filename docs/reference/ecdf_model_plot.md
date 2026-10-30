# ECDF Plot Comparing Model Strategies (fixed prior)

ECDF (+ optional density) comparing sensitivity model strategies – All
studies, Excluding High RoB, PET-PEESE, Mixture Model, BMA – for a
single prior. Colour represents the model strategy/section.

## Usage

``` r
ecdf_model_plot(
  model,
  data,
  priors,
  measure,
  prior_to_plot = NULL,
  study_var = NULL,
  rob_var = NULL,
  exclude_high_rob = FALSE,
  incl_common_effect = FALSE,
  incl_random_effect = TRUE,
  incl_pet_peese = FALSE,
  incl_selection_weight = FALSE,
  incl_robma = FALSE,
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

- priors:

  Named list of prior specs; each element a list with `name`,
  `mu_prior`, `tau_prior`.

- measure:

  Effect measure ("OR","RR","HR","IRR","MD","SMD").

- prior_to_plot:

  Name (in `priors`) of the prior to display. Defaults to the first
  prior.

- study_var, rob_var:

  Optional unquoted column names (study id / RoB).

- exclude_high_rob, incl_common_effect, incl_random_effect,
  incl_pet_peese, incl_selection_weight, incl_robma:

  Logical section toggles.

- robma_fit:

  A fitted RoBMA object (required if `incl_robma`).

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
