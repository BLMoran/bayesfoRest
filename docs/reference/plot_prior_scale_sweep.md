# Plot a prior-scale sweep (mu vs prior SD)

Plot a prior-scale sweep (mu vs prior SD)

## Usage

``` r
plot_prior_scale_sweep(
  sweep,
  prior_mean = 0,
  x_label = "mu-prior SD (log scale)",
  y_label = "Pooled effect",
  font = NULL
)
```

## Arguments

- sweep:

  Output of
  [`prior_scale_sweep()`](https://blmoran.github.io/bayesfoRest/reference/prior_scale_sweep.md).

- prior_mean:

  The prior mean used (drawn as a reference line).

- x_label, y_label:

  Axis labels.

- font:

  Font family.

## Value

A ggplot.
