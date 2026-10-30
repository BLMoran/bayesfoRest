# Extract Pooled Effect Summary

Extracts the pooled (overall) effect estimate and credible interval from
a fitted model. Values are on the model scale (log scale for ratio
measures).

## Usage

``` r
extract_pooled_summary(model, prob = 0.95, ...)
```

## Arguments

- model:

  A fitted model object.

- prob:

  Probability mass for the credible interval (default 0.95).

- ...:

  Additional arguments passed to methods.

## Value

A named numeric vector: `estimate`, `lower`, `upper`.
