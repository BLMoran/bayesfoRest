# Plot a tau-prior scale sweep (mu and tau vs tau-prior SD)

Plot a tau-prior scale sweep (mu and tau vs tau-prior SD)

## Usage

``` r
plot_tau_scale_sweep(
  sweep,
  mu_label = "Pooled effect (μ)",
  tau_label = "Heterogeneity (τ)",
  null_value = 0,
  font = NULL
)
```

## Arguments

- sweep:

  Output of
  [`tau_scale_sweep()`](https://blmoran.github.io/bayesfoRest/reference/tau_scale_sweep.md).

- mu_label, tau_label:

  Facet labels for the two panels.

- null_value:

  Reference line for mu (default 0).

- font:

  Font family.

## Value

A ggplot (two stacked panels: mu and tau vs tau-prior SD).
