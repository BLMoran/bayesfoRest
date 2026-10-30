# Generate JAGS Model Files for Bayesian Meta-Analysis

Writes a JAGS model file for common meta-analysis models (binary,
continuous, count outcomes) with a user-specified prior on the
between-study standard deviation (sigma). These model files can be
passed directly to
[`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html).

## Usage

``` r
jags_model_file(
  outcome = c("bin", "ctns", "count"),
  model,
  type = c("ran", "fix"),
  sigma_prior = jags_prior_halfnormal(0.5),
  path = NULL
)

jags_prior_halfnormal(scale = 0.5)

jags_prior_halfcauchy(scale = 0.5)

jags_prior_uniform(upper = 5)
```

## Arguments

- outcome:

  Character: `"bin"`, `"ctns"`, or `"count"`.

- model:

  Character: model specification string.

  - Binary: `"std.norm"`, `"std.dt"`, `"reg.norm"`, `"reg.dt"`

  - Continuous: `"std.ta"`, `"std.mv"`, `"reg.ta"`, `"reg.mv"`

  - Count: `"std"`, `"std.unif"`, `"std.hc"`, `"reg"`, `"reg.unif"`,
    `"reg.hc"`

- type:

  Character: `"fix"` (common-effects) or `"ran"` (random-effects).

- sigma_prior:

  Character: a block of JAGS code defining `sigma` and `tau`. Must
  define both variables. Defaults to a half-normal(0, 0.5) prior. Use
  the helper functions:

  - `jags_prior_halfnormal(scale)` (default)

  - `jags_prior_halfcauchy(scale)`

  - `jags_prior_uniform(upper)`

- path:

  File path to write the model to. Defaults to a tempfile.

- scale:

  Scale parameter of the half-Cauchy. Default 0.5.

- upper:

  Upper bound of the uniform. Default 5.

## Value

The file path (invisibly), for passing to
[`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html).

## Functions

- `jags_prior_halfnormal()`: Half-normal prior on sigma:
  `sigma ~ N(0, scale^2) T(0,)`

- `jags_prior_halfcauchy()`: Half-Cauchy prior on sigma via the
  scale-mixture representation.

- `jags_prior_uniform()`: Uniform prior on sigma.
