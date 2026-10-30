# Generate Sensitivity Analysis Plot for Bayesian Meta-Analysis

Creates a comprehensive sensitivity analysis visualization showing how
meta-analytic estimates vary across different prior specifications,
analysis approaches, and publication bias adjustments. Supports all
model backends (brms, MetaStan, bayesmeta, rstan, cmdstanr) and
integrates RoBMA for bias-adjusted sensitivity analyses.

## Usage

``` r
sensitivity_plot(
  model,
  data,
  priors,
  measure,
  prior_order = NULL,
  model_order = NULL,
  rob_var = NULL,
  exclude_high_rob = FALSE,
  incl_common_effect = FALSE,
  incl_random_effect = TRUE,
  incl_selection_weight = FALSE,
  incl_pet_peese = FALSE,
  incl_robma = FALSE,
  robma_fit = NULL,
  seed = TRUE,
  add_probs = FALSE,
  null_value = NULL,
  null_range = NULL,
  add_null_range = FALSE,
  color_null_range = "#77bb41",
  label_control = "Control",
  label_intervention = "Intervention",
  title = NULL,
  subtitle = NULL,
  title_align = "left",
  xlim = NULL,
  x_breaks = NULL,
  color_palette = NULL,
  color_overall_posterior = "dodgerblue",
  color_overall_posterior_outline = "blue",
  split_color_by_null = FALSE,
  color_favours_control = "firebrick",
  color_favours_intervention = "dodgerblue",
  plot_width = 4,
  font = NULL
)
```

## Arguments

- model:

  A fitted model object (brmsfit, meta_stan, bayesmeta, stanreg,
  stanfit, or CmdStanMCMC).

- data:

  A data frame containing the study data.

- priors:

  A named list of prior specifications. Each element should be a list
  with `name` (display label), `mu_prior` (a `bf_prior` object), and
  optionally `tau_prior` (a `bf_prior` object). See
  [`bf_normal`](https://blmoran.github.io/bayesfoRest/reference/bf_priors.md)
  for constructors.

- measure:

  Character. Effect measure: "OR", "RR", "HR", "IRR", "MD", or "SMD".

- prior_order:

  Optional character vector of prior keys specifying display order.
  Defaults to the order in `priors`.

- model_order:

  Optional character vector of section labels specifying display order.
  Overrides the default ordering.

- rob_var:

  Optional. Name of risk of bias variable (unquoted).

- exclude_high_rob:

  Logical. Include sensitivity excluding high RoB studies. Default
  FALSE.

- incl_common_effect:

  Logical. Include common-effect model section. Default FALSE. Only
  applies when the base model is random-effects.

- incl_random_effect:

  Logical. Include random-effects model section. Default TRUE.

- incl_selection_weight:

  Logical. Include RoBMA selection (weight) model rows. Default FALSE.
  Requires `robma_fit`.

- incl_pet_peese:

  Logical. Include RoBMA PET-PEESE rows. Default FALSE. Requires
  `robma_fit`.

- incl_robma:

  Logical. Include full RoBMA sensitivity rows (conditional and
  model-averaged). Default FALSE. Requires `robma_fit`.

- robma_fit:

  A single fitted RoBMA object (from
  [`RoBMA::RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.html))
  or a named list of RoBMA objects (one per prior level). Required when
  any RoBMA option is TRUE.

- seed:

  Integer seed for reproducibility. Default TRUE (uses model seed).

- add_probs:

  Logical. Add probability columns. Default FALSE.

- null_value:

  Numeric. Null value for hypothesis testing.

- null_range:

  Numeric vector of length 2 or single value. ROPE.

- add_null_range:

  Logical. Show null range region. Default FALSE.

- color_null_range:

  Character. Color for null range. Default "#77bb41".

- label_control:

  Character. Control group label. Default "Control".

- label_intervention:

  Character. Intervention label. Default "Intervention".

- title:

  Optional title.

- subtitle:

  Optional subtitle.

- title_align:

  Character. "left", "center", or "right".

- xlim:

  Optional x-axis limits.

- x_breaks:

  Optional x-axis breaks.

- color_palette:

  Optional named vector of colors for sections.

- color_overall_posterior:

  Fill color. Default "dodgerblue".

- color_overall_posterior_outline:

  Outline color. Default "blue".

- split_color_by_null:

  Logical. Color by null. Default FALSE.

- color_favours_control:

  Character. Default "firebrick".

- color_favours_intervention:

  Character. Default "dodgerblue".

- plot_width:

  Numeric. Width ratio for density plot. Default 4.

- font:

  Optional font family.

## Value

A patchwork object combining tables and density plots.

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

sensitivity_plot(
  model   = fit,
  data    = dat,
  priors  = priors,
  measure = "OR",
  incl_common_effect = TRUE,
  exclude_high_rob   = TRUE,
  rob_var            = rob_rating
)
} # }
```
