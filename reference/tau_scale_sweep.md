# Continuous tau-prior scale sweep

Continuous tau-prior scale sweep

## Usage

``` r
tau_scale_sweep(
  fit,
  scales = c(0.25, 0.5, 0.75, 1, 1.5, 3),
  mu_par = NULL,
  refresh = 0
)
```

## Arguments

- fit:

  A brms / rstanarm meta-analysis fit.

- scales:

  Numeric vector of half-normal scales for the tau prior.

- mu_par:

  Pooled-effect draw name (auto-detected).

- refresh:

  Passed to refits.

## Value

A tibble: tau_prior_sd, mu_median, mu_lo, mu_hi, pr_benefit, tau_median,
tau_lo, tau_hi, tau_prior_pctile, conflict.
