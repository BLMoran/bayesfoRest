# Using bayesfoRest with brms

brms is the most fully-featured backend in bayesfoRest. Every plotting
function works with a `brmsfit`, including the ones that need to refit
the model: subgroup analyses and prior sensitivity
([`sensitivity_plot()`](https://blmoran.github.io/bayesfoRest/reference/sensitivity_plot.md),
[`ecdf_model_plot()`](https://blmoran.github.io/bayesfoRest/reference/ecdf_model_plot.md),
[`ecdf_prior_plot()`](https://blmoran.github.io/bayesfoRest/reference/ecdf_prior_plot.md)).
If you’re not sure which backend to start with, start with brms.

## Fit the model

bayesfoRest expects a two-column effect-and-standard-error outcome, fit
with [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).
Prepare the effect sizes with
[`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
first:

``` r

library(metafor)
library(brms)
library(bayesfoRest)

binary_es <- escalc(
  measure = "OR",
  ai = Event_Intervention, bi = Outcome_Intervention_No,
  ci = Event_Control, di = Outcome_Control_No,
  n1i = N_Intervention, n2i = N_Control,
  data = binary_outcome) |>
  dplyr::mutate(sei = sqrt(vi))
```

A random-effect (hierarchical) model includes `(1 | Author)`; a
common-effect model drops it:

``` r

# Random-effect
model_re <- brm(
  yi | se(sei) ~ 1 + (1 | Author),
  family = gaussian(),
  data = binary_es,
  prior = prior(normal(0, 1), class = Intercept) +
          prior(cauchy(0, 0.5), class = sd),
  iter = 4000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = 0.99),
  silent = 2)

# Common-effect
model_ce <- brm(
  yi | se(sei) ~ 1,
  family = gaussian(),
  data = binary_es,
  prior = prior(normal(0, 1), class = Intercept),
  iter = 4000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = 0.99),
  silent = 2)
```

[`bayes_forest()`](https://blmoran.github.io/bayesfoRest/reference/bayes_forest.md)
(and every other bayesfoRest function) auto-detects the pooling type
from the model by checking for a `sd_` draw – you don’t need to pass
`model_type` unless you want to override it. See
[`vignette("bayesfoRest")`](https://blmoran.github.io/bayesfoRest/articles/bayesfoRest.md)
for the full forest plot walkthrough; this vignette focuses on what’s
brms-specific.

## Subgroup analysis

With `subgroup = TRUE`, bayesfoRest needs a per-subgroup pooled
estimate. For brms it gets one automatically:
`refit_for_subgroup.brmsfit()` calls
`stats::update(model, newdata = <subgroup rows>, recompile = FALSE)` for
each subgroup, reusing the compiled Stan model so this stays fast.

``` r

bayes_forest(
  model = model_re,
  data = binary_es,
  measure = "OR",
  studyvar = Author,
  yearvar = Year,
  n_ctrl = N_Control,
  n_int = N_Intervention,
  event_ctrl = Event_Control,
  event_int = Event_Intervention,
  subgroup = TRUE,
  subgroupvar = Subgroup,
  xlim = c(0.1, 3.5)
)
```

Backends without a `refit_for_subgroup()` method fall back to
approximating subgroup effects from the overall model, with a warning.
brms (along with bayesmeta, MetaStan, rstanarm and bayesma) refits
properly.

## Overall and funnel plots

``` r

overall_plot(
  model = model_re,
  data = binary_es,
  measure = "OR",
  add_tau = TRUE,
  add_pred = TRUE,
  incl_mu_prior = TRUE,
  incl_tau_prior = TRUE
)
```

`incl_mu_prior`/`incl_tau_prior` read the priors directly from the brms
fit (via
[`brms::prior_summary()`](https://mc-stan.org/rstantools/reference/prior_summary.html))
and draw them on the same density scale as the posteriors – no need to
re-specify them. This is currently only automatic for brms and
bayesmeta; other backends need `mu_prior`/`tau_prior` passed in
explicitly (e.g. `mu_prior = bf_normal(0, 1)`).

``` r

funnel_plot(
  model = model_re,
  data = binary_es,
  measure = "OR",
  studyvar = Author,
  yearvar = Year,
  yi_col = "yi",
  sei_col = "sei"
)
```

## Prior sensitivity

brms can fit every prior family `bf_*()` offers –
[`bf_normal()`](https://blmoran.github.io/bayesfoRest/reference/bf_priors.md),
[`bf_student_t()`](https://blmoran.github.io/bayesfoRest/reference/bf_priors.md),
[`bf_half_normal()`](https://blmoran.github.io/bayesfoRest/reference/bf_priors.md),
[`bf_half_cauchy()`](https://blmoran.github.io/bayesfoRest/reference/bf_priors.md),
[`bf_half_t()`](https://blmoran.github.io/bayesfoRest/reference/bf_priors.md),
[`bf_exponential()`](https://blmoran.github.io/bayesfoRest/reference/bf_priors.md),
[`bf_lognormal()`](https://blmoran.github.io/bayesfoRest/reference/bf_priors.md),
[`bf_uniform()`](https://blmoran.github.io/bayesfoRest/reference/bf_priors.md),
and the empirical
[`bf_turner()`](https://blmoran.github.io/bayesfoRest/reference/bf_turner.md)/[`bf_rhodes()`](https://blmoran.github.io/bayesfoRest/reference/bf_turner.md)
priors – so
[`sensitivity_plot()`](https://blmoran.github.io/bayesfoRest/reference/sensitivity_plot.md),
[`ecdf_model_plot()`](https://blmoran.github.io/bayesfoRest/reference/ecdf_model_plot.md)
and
[`ecdf_prior_plot()`](https://blmoran.github.io/bayesfoRest/reference/ecdf_prior_plot.md)
all work without restriction. Write the priors you want to compare with
the `bf_*()` helpers:

``` r

priors <- list(
  vague = list(
    name      = "Vague",
    mu_prior  = bf_normal(0, 10),
    tau_prior = bf_half_cauchy(0, 1)),
  weakreg = list(
    name      = "Weakly Regularising",
    mu_prior  = bf_normal(0, 1),
    tau_prior = bf_half_cauchy(0, 0.5)),
  informative = list(
    name      = "Informative",
    mu_prior  = bf_normal(0, 0.5),
    tau_prior = bf_half_cauchy(0, 0.3))
)
```

[`sensitivity_plot()`](https://blmoran.github.io/bayesfoRest/reference/sensitivity_plot.md)
refits `model_re` under each prior and draws the pooled effect from
every section you ask for. `incl_random_effect` (default `TRUE`) keeps
the original random-effect sections; `incl_common_effect = TRUE` adds a
common-effect refit for the same priors; `exclude_high_rob = TRUE` adds
a section excluding studies with a high overall risk of bias (needs
`rob_var`):

``` r

sensitivity_plot(
  model = model_re,
  data = binary_es,
  priors = priors,
  measure = "OR",
  rob_var = Overall,
  incl_random_effect = TRUE,
  incl_common_effect = TRUE,
  exclude_high_rob = TRUE,
  xlim = c(0.25, 1.5)
)
```

[`ecdf_model_plot()`](https://blmoran.github.io/bayesfoRest/reference/ecdf_model_plot.md)
fixes one prior and compares strategies (sections);
[`ecdf_prior_plot()`](https://blmoran.github.io/bayesfoRest/reference/ecdf_prior_plot.md)
fixes one section and compares priors:

``` r

ecdf_model_plot(
  model = model_re,
  data = binary_es,
  priors = priors["weakreg"],
  measure = "OR",
  rob_var = Overall,
  incl_common_effect = TRUE,
  exclude_high_rob = TRUE,
  xlim = c(0.25, 1.5)
)

ecdf_prior_plot(
  model = model_re,
  data = binary_es,
  priors = priors,
  measure = "OR",
  xlim = c(0.25, 1.5)
)
```

See `vignette("sensitivity-analysis")` for publication-bias-adjusted
sensitivity sections built from a RoBMA fit (`incl_robma`,
`incl_pet_peese`, `incl_selection_weight`).
