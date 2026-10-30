# Prior Distribution Constructors for Sensitivity Analysis

Backend-agnostic prior specification functions for use with
[`sensitivity_plot`](https://blmoran.github.io/bayesfoRest/reference/sensitivity_plot.md),
[`ecdf_model_plot`](https://blmoran.github.io/bayesfoRest/reference/ecdf_model_plot.md),
[`ecdf_prior_plot`](https://blmoran.github.io/bayesfoRest/reference/ecdf_prior_plot.md)
and
[`overall_plot`](https://blmoran.github.io/bayesfoRest/reference/overall_plot.md).
These create lightweight prior objects that are translated to the format
each model backend expects. A prior a backend cannot fit exactly is an
error, never a silent substitution.

## Usage

``` r
bf_normal(mean = 0, sd = 1)

bf_student_t(df = 3, location = 0, scale = 1)

bf_half_normal(mean = 0, sd = 1)

bf_half_cauchy(location = 0, scale = 1)

bf_half_t(df = 3, location = 0, scale = 1)

bf_exponential(rate = 1)

bf_lognormal(meanlog = 0, sdlog = 1)

bf_uniform(lower = 0, upper = 1)
```

## Arguments

- mean, location:

  Prior mean or location parameter.

- sd, scale:

  Prior standard deviation or scale parameter.

- df:

  Degrees of freedom.

- rate:

  Rate parameter (for `bf_exponential` only).

- meanlog, sdlog:

  Mean and standard deviation of \\\log(\tau)\\ (for `bf_lognormal`
  only).

- lower, upper:

  Bounds of the uniform prior. For \\\tau\\, `lower` must be 0.

## Value

An object of class `"bf_prior"` containing the distribution family and
parameters.

## Details

Priors for the pooled effect (\\\mu\\) are on the model scale (the log
scale for OR, RR, HR and IRR). Priors for the heterogeneity (\\\tau\\)
are on the standard deviation of the study effects.

**Priors for \\\mu\\**

- `bf_normal(mean, sd)`:

  Normal distribution.

- `bf_student_t(df, location, scale)`:

  Student-t distribution: a heavier-tailed, weakly informative
  alternative to the normal. `df = 1` is a Cauchy distribution.

**Priors for \\\tau\\**

- `bf_half_normal(mean, sd)`:

  Half-normal distribution (truncated at 0).

- `bf_half_cauchy(location, scale)`:

  Half-Cauchy distribution.

- `bf_half_t(df, location, scale)`:

  Half-Student-t distribution.

- `bf_exponential(rate)`:

  Exponential distribution.

- `bf_lognormal(meanlog, sdlog)`:

  Log-normal distribution. Use
  [`bf_turner`](https://blmoran.github.io/bayesfoRest/reference/bf_turner.md)
  for the empirical priors of Turner et al.

- `bf_uniform(0, upper)`:

  Uniform distribution from 0 to `upper`.

- [`bf_rhodes()`](https://blmoran.github.io/bayesfoRest/reference/bf_turner.md):

  Empirical log-t priors of Rhodes et al.

**Backend support**

Sensitivity refits stop with an error if the model's backend cannot fit
a prior exactly:

|                             |      |           |                  |
|-----------------------------|------|-----------|------------------|
| Prior                       | brms | bayesmeta | MetaStan         |
| `bf_normal` (\\\mu\\)       | yes  | yes       | yes              |
| `bf_student_t` (\\\mu\\)    | yes  | no        | no               |
| `bf_half_normal`            | yes  | yes       | yes (mean 0)     |
| `bf_half_cauchy`            | yes  | yes       | yes (location 0) |
| `bf_half_t`                 | yes  | yes       | no               |
| `bf_exponential`            | yes  | yes       | no               |
| `bf_lognormal`, `bf_turner` | yes  | yes       | no               |
| `bf_uniform`                | yes  | yes       | yes              |
| `bf_rhodes`                 | yes  | yes       | no               |

Priors cannot be changed for rstan, cmdstanr, rstanarm, bayesma, RoBMA
or JAGS models.

## Examples

``` r
if (FALSE) { # \dontrun{
priors <- list(
  vague = list(
    name      = "Vague",
    mu_prior  = bf_normal(0, 10),
    tau_prior = bf_half_cauchy(0, 1)
  ),
  weakreg = list(
    name      = "Weakly Regularising",
    mu_prior  = bf_student_t(3, 0, 1),
    tau_prior = bf_half_normal(0, 0.5)
  ),
  empirical = list(
    name      = "Empirical (Turner)",
    mu_prior  = bf_normal(0, 1),
    tau_prior = bf_turner("all-cause mortality",
                          "pharmacological", "placebo / control")
  )
)

sensitivity_plot(model, data, priors = priors, measure = "OR")
} # }
```
