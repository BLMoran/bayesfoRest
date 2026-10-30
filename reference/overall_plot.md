# Create posterior plots for Bayesian meta-analysis

Produces density plots of the pooled effect (\\\mu\\), heterogeneity
(\\\tau\\), and/or prediction distribution from a fitted model. Supports
all bayesfoRest backends: brms, MetaStan, bayesmeta, rstan, cmdstanr,
and rstanarm.

## Usage

``` r
overall_plot(
  model,
  data,
  measure,
  add_tau = FALSE,
  add_pred = FALSE,
  add_summary = TRUE,
  incl_mu_prior = FALSE,
  incl_tau_prior = FALSE,
  mu_prior = NULL,
  tau_prior = NULL,
  color_prior = "grey40",
  color_pred_posterior = "orange",
  color_pred_posterior_outline = "darkorange2",
  null_value = NULL,
  null_range = NULL,
  add_null_range = FALSE,
  color_null_range = "#77bb41",
  label_control = "Control",
  label_intervention = "Intervention",
  title = NULL,
  subtitle = NULL,
  title_align = "left",
  x_label = NULL,
  mu_xlim = NULL,
  tau_xlim = NULL,
  pred_xlim = NULL,
  x_breaks = NULL,
  tau_breaks = NULL,
  pred_breaks = NULL,
  color_overall_posterior = "dodgerblue",
  color_overall_posterior_outline = "blue",
  split_color_by_null = FALSE,
  color_favours_control = "firebrick",
  color_favours_intervention = "dodgerblue",
  tau_posterior_color = "#77bb41",
  tau_posterior_outline = NULL,
  tau_slab_scale = 0.7,
  font = NULL,
  plot_arrangement = "vertical",
  ...
)
```

## Arguments

- model:

  A fitted model object (brmsfit, meta_stan, bayesmeta, stanreg,
  stanfit, or CmdStanMCMC).

- data:

  Data frame with study information. Must include an `Author` column and
  effect size / SE columns appropriate for the backend.

- measure:

  Character string for effect measure ("OR", "RR", "HR", "IRR", "MD",
  "SMD").

- add_tau:

  Logical, whether to include the heterogeneity (tau) plot.

- add_pred:

  Logical, whether to include the prediction distribution plot.

- add_summary:

  Logical, whether to annotate the mu plot with a summary of the effect
  estimate, 95\\ benefit/harm. When `null_range` is provided, also shows
  Pr(Benefit \> \\\delta\\) and Pr(Harm \> \\\delta\\). Default `TRUE`.

- incl_mu_prior:

  Logical, whether to draw the prior for \\\mu\\ as a dashed line behind
  the posterior. Default `FALSE`.

- incl_tau_prior:

  Logical, whether to draw the prior for \\\tau\\ as a dashed line
  behind the posterior (needs `add_tau = TRUE`). Default `FALSE`.

- mu_prior, tau_prior:

  Optional priors to draw, as `bf_prior` objects (e.g.
  `bf_normal(0, 1)`, `bf_half_cauchy(0, 0.5)`) or density functions on
  the model scale (log scale for ratio measures). When `NULL`, the
  priors are read from brms and bayesmeta models. Other backends do not
  store their priors, so supply them here.

- color_prior:

  Color for the prior lines.

- color_pred_posterior:

  Color for the prediction posterior distribution.

- color_pred_posterior_outline:

  Outline color for the prediction posterior distribution.

- null_value:

  Numeric, the null value for the effect (default from
  `get_measure_properties`).

- null_range:

  Range for the null region: a numeric vector of length 2 (e.g.,
  `c(0.9, 1.1)`) or a single half-width around `null_value`. When `NULL`
  and `add_null_range = TRUE`, defaults to 0.9 to 1.1 for ratio measures
  and -0.1 to 0.1 for SMD.

- add_null_range:

  Logical, whether to add shaded null range region.

- color_null_range:

  Color for the null range shading.

- label_control:

  Character string for control group label.

- label_intervention:

  Character string for intervention group label.

- title:

  Plot title.

- subtitle:

  Plot subtitle.

- title_align:

  Alignment for title ("left", "center", "right").

- x_label:

  Character string for the mu/prediction plot x-axis label.

- mu_xlim:

  Numeric vector of length 2 for mu plot x-axis limits.

- tau_xlim:

  Numeric vector of length 2 for tau plot x-axis limits.

- pred_xlim:

  Numeric vector of length 2 for prediction plot x-axis limits.

- x_breaks:

  Numeric vector for x-axis breaks on mu plot.

- tau_breaks:

  Numeric vector for x-axis breaks on tau plot.

- pred_breaks:

  Numeric vector for x-axis breaks on prediction plot.

- color_overall_posterior:

  Color for the overall posterior distribution.

- color_overall_posterior_outline:

  Outline color for the posterior distribution.

- split_color_by_null:

  Logical, whether to split posterior color by null value.

- color_favours_control:

  Color for region favoring control.

- color_favours_intervention:

  Color for region favoring intervention.

- tau_posterior_color:

  Base color for tau posterior.

- tau_posterior_outline:

  Outline color for tau posterior.

- tau_slab_scale:

  Numeric, scaling factor for tau slab height relative to mu (default
  0.7).

- font:

  Font family for text elements.

- plot_arrangement:

  Character, "vertical" (stacked) or "horizontal" (side by side).

- ...:

  Additional arguments passed to `extract_forest_draws` (e.g.,
  `theta_par`, `tau_par`, `study_par` for stanfit/CmdStanMCMC models).

## Value

A ggplot object (or combined plot if `add_tau = TRUE` or
`add_pred = TRUE`).

## Details

Priors and posteriors are drawn on the same density scale, so a vague
prior appears almost flat and a prior that dominates the posterior
appears taller than it.
