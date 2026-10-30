# Assess prior dominance in a Bayesian meta-analysis (backend-agnostic)

Assess prior dominance in a Bayesian meta-analysis (backend-agnostic)

## Usage

``` r
prior_dominance(
  fit,
  prior_draws = NULL,
  mu_par = NULL,
  tau_par = NULL,
  tau_grid = c(0.5, 0.75, 1.5),
  refresh = 0
)
```

## Arguments

- fit:

  A fitted model. Posterior draws are read via
  [`posterior::as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html),
  so any backend that supports it works.

- prior_draws:

  Optional PRIOR draws (data frame / draws object) with the mu and tau
  columns. If supplied it is used directly (portable path for cmdstanr /
  rstan / JAGS). If NULL, priors are sampled from the fit.

- mu_par, tau_par:

  Column names for the pooled effect and between-study SD. Auto-detected
  when NULL.

- tau_grid:

  Half-normal scales for the tau-prior sensitivity refits (brms /
  rstanarm only). NULL to skip.

- refresh:

  Passed to refits.

## Value

An object of class `prior_dominance`.
