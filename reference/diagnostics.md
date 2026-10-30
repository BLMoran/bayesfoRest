# Single-Page MCMC Diagnostics

Single-Page MCMC Diagnostics

## Usage

``` r
diagnostics(
  object,
  pars = NULL,
  labels = NULL,
  study_labels = NULL,
  model_label = NULL,
  ndraws = 100,
  font = NULL,
  ...
)
```

## Arguments

- object:

  A fitted model: brmsfit, stanfit (rstan / MetaStan), CmdStanMCMC
  (cmdstanr), or stanreg.

- pars:

  Character vector of parameters for the trace / autocorrelation panels
  (given as ORIGINAL parameter names). Default selects the
  population-level effect and between-study SD.

- labels:

  Optional named vector mapping original parameter names to display
  labels, e.g.
  `c(b_Intercept = "\u03bc", sd_Author__Intercept = "\u03c4")`.
  Overrides the automatic mu / tau defaults.

- study_labels:

  Optional named vector mapping study keys to display labels, e.g. as
  built by
  [`make_study_labels()`](https://blmoran.github.io/bayesfoRest/reference/make_study_labels.md).

- model_label:

  Optional character string naming the model, shown in the diagnostics
  page title.

- ndraws:

  Draws used in the posterior predictive panel.

- font:

  Font family for panel text.

- ...:

  Unused.

## Value

A patchwork object (the assembled diagnostic page).
