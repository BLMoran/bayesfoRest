# Detect Model Backend and Pooling Type

Inspects a fitted model object and returns its backend type and whether
it uses common-effects or random-effects pooling.

## Usage

``` r
detect_model_info(model, ...)
```

## Arguments

- model:

  A fitted model object.

- ...:

  Additional arguments passed to methods.

## Value

A list with elements:

- backend:

  Character string identifying the backend.

- pooling:

  Character string: "random_effect" or "common_effect".
