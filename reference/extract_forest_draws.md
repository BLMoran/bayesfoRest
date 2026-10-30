# Extract Posterior Draws in Canonical Format

Extracts study-level and pooled posterior draws from a fitted model,
returning them in the canonical format expected by bayesfoRest plotting
and table functions.

## Usage

``` r
extract_forest_draws(model, data, ...)
```

## Arguments

- model:

  A fitted model object.

- data:

  Data frame with study information. Must include an `Author` column.

- ...:

  Additional arguments passed to methods (e.g., `theta_par`, `tau_par`,
  `study_par` for stanfit/CmdStanMCMC models).

## Value

A tibble with columns: `.draw`, `Author`, `theta_study`, `eta`, `tau`.
