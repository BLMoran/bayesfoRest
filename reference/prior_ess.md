# Prior effective sample size for the pooled effect (study-equivalents)

Prior effective sample size for the pooled effect (study-equivalents)

## Usage

``` r
prior_ess(
  fit,
  data = NULL,
  typical_se = NULL,
  prior_draws = NULL,
  mu_par = NULL
)
```

## Arguments

- fit:

  A fitted meta-analysis model.

- data:

  The study-level data (needs the sei column) to define a typical study
  precision. Alternatively pass `typical_se` directly.

- typical_se:

  Optional numeric: the representative within-study SE. If NULL,
  computed as the pooled (harmonic-mean-precision) SE from `data$sei`.

- prior_draws, mu_par:

  As in
  [`prior_dominance()`](https://blmoran.github.io/bayesfoRest/reference/prior_dominance.md).

## Value

A list: prior_sd, typical_se, prior_ess (study-equivalents), k (actual
studies), and the prior's share of total precision.
