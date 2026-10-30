# Continuous mu-prior scale sweep

Continuous mu-prior scale sweep

## Usage

``` r
prior_scale_sweep(
  fit,
  scales = c(0.1, 0.25, 0.5, 1, 2, 5),
  prior_mean = 0,
  mu_par = NULL,
  refresh = 0
)
```

## Arguments

- fit:

  A brms / rstanarm meta-analysis fit.

- scales:

  Numeric vector of mu-prior SDs to sweep.

- prior_mean:

  Mean of the mu prior held fixed across the sweep (default 0). Set to
  your informative-prior mean to test how hard a tight informative prior
  pulls the estimate.

- mu_par:

  Pooled-effect draw name (auto-detected).

- refresh:

  Passed to refits.

## Value

A tibble: prior_sd, mu_median, mu_lo, mu_hi, pr_benefit.
