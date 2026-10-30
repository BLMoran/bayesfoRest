# ECDF Plot Comparing Priors (fixed section)

ECDF (+ optional density) comparing the posterior across priors for a
single analysis section (default "All studies"). Colour represents the
prior. This is the prior-sensitivity companion to
[`ecdf_model_plot`](https://blmoran.github.io/bayesfoRest/reference/ecdf_model_plot.md).

## Usage

``` r
ecdf_prior_plot(
  model,
  data,
  priors,
  measure,
  section_to_plot = "Random Effects",
  prior_order = NULL,
  study_var = NULL,
  rob_var = NULL,
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
  font = NULL,
  ...
)
```

## Arguments

- model:

  A fitted `brmsfit` meta-analysis model.

- data:

  Study-level data used to fit `model`.

- priors:

  Named list of prior specs; each element a list with `name`,
  `mu_prior`, `tau_prior`. The display order of the legend/colours
  follows the order of this list (or `prior_order`).

- measure:

  Effect measure ("OR","RR","HR","IRR","MD","SMD").

- section_to_plot:

  Which analysis section to display. One of "All studies", "Excluding
  High RoB", "PET-PEESE", "Mixture Model", "BMA". Default "All studies".

- prior_order:

  Optional character vector giving the order of priors in the legend.
  Defaults to `names(priors)`.

- study_var, rob_var:

  Optional unquoted column names (study id / RoB), needed only if the
  chosen `section_to_plot` requires them.

- prob_reference:

  "null" or "null_range".

- null_value, null_range, add_null_range, color_null_range:

  Null / ROPE options.

- label_control, label_intervention:

  Arm labels for the density panel.

- title, subtitle, x_label, xlim, x_breaks, color_palette, show_density,
  font:

  Display options.

- ...:

  Additional section toggles (e.g. `exclude_high_rob`, `incl_pet_peese`,
  `incl_selection_weight`, `incl_robma`) passed on to
  `.build_ecdf_draws()`. Set the one corresponding to `section_to_plot`
  (other than "All studies") to TRUE so those draws are built.

## Value

A ggplot (or patchwork) object.

## See also

[`ecdf_model_plot`](https://blmoran.github.io/bayesfoRest/reference/ecdf_model_plot.md)
for comparing model strategies at a fixed prior.
