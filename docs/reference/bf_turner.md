# Empirical Heterogeneity Priors

Informative priors for \\\tau\\ from large analyses of Cochrane reviews.
They are most useful when there are few studies, where \\\tau\\ is
poorly identified by the data.

`bf_turner()` gives the log-normal priors of Turner et al. (2015) for
log odds ratios (binary outcomes). `bf_rhodes()` gives the log-t priors
of Rhodes et al. (2015) for standardised mean differences (continuous
outcomes). The parameters come from
[`bayesmeta::TurnerEtAlPrior()`](https://rdrr.io/pkg/bayesmeta/man/TurnerEtAlPrior.html)
and
[`bayesmeta::RhodesEtAlPrior()`](https://rdrr.io/pkg/bayesmeta/man/RhodesEtAlPrior.html),
so the bayesmeta package must be installed. See those functions for the
outcome and comparison categories.

## Usage

``` r
bf_turner(
  outcome = NA,
  comparator1 = "pharmacological",
  comparator2 = "placebo / control"
)

bf_rhodes(
  outcome = NA,
  comparator1 = "pharmacological",
  comparator2 = "placebo / control",
  area = "other"
)
```

## Arguments

- outcome:

  Outcome type, e.g. `"all-cause mortality"` for `bf_turner()` or
  `"mental health outcome"` for `bf_rhodes()`. `NA` gives the prior for
  a general outcome.

- comparator1, comparator2:

  Intervention types compared: any of `"pharmacological"`,
  `"non-pharmacological"` and `"placebo / control"`.

- area:

  Medical area for `bf_rhodes()`: `"other"`, `"respiratory"` or
  `"cancer"`.

## Value

A `bf_prior` for \\\tau\\. See
[`bf_priors`](https://blmoran.github.io/bayesfoRest/reference/bf_priors.md)
for which backends support it.

## References

Turner RM, Jackson D, Wei Y, Thompson SG, Higgins JPT (2015). Predictive
distributions for between-study heterogeneity and simple methods for
their application in Bayesian meta-analysis. *Statistics in Medicine*,
34(6), 984–998.

Rhodes KM, Turner RM, Higgins JPT (2015). Predictive distributions were
developed for the extent of heterogeneity in meta-analyses of continuous
outcome data. *Journal of Clinical Epidemiology*, 68(1), 52–60.

## Examples

``` r
if (FALSE) { # \dontrun{
bf_turner("all-cause mortality", "pharmacological", "placebo / control")
bf_rhodes("mental health outcome", "non-pharmacological", "placebo / control")
} # }
```
