# Create ECDF Plot for Sensitivity Analysis

Creates an empirical cumulative distribution function (ECDF) plot
showing posterior distributions across different sensitivity analyses
(All studies, Excluding High RoB, PET-PEESE, Mixture Model, BMA).

## Usage

``` r
ecdf_plot(
  model,
  data,
  priors,
  measure,
  study_var = NULL,
  rob_var = NULL,
  exclude_high_rob = FALSE,
  incl_pet_peese = FALSE,
  pet_peese_direction = "negative",
  pet_peese_threshold = 0.1,
  incl_mixture = FALSE,
  incl_bma = FALSE,
  model_bma = NULL,
  prob_reference = "null",
  null_value = NULL,
  null_range = NULL,
  add_null_range = FALSE,
  color_null_range = "#77bb41",
  label_control = "Control",
  label_intervention = "Intervention",
  title = NULL,
  subtitle = NULL,
  xlim = NULL,
  x_breaks = NULL,
  color_palette = NULL,
  show_density = TRUE,
  prior_to_plot = "weakreg",
  font = NULL
)
```

## Arguments

- model:

  A fitted brmsfit object from a Bayesian meta-analysis

- data:

  A data frame containing the study data used to fit the model

- priors:

  A named list of brms prior specifications (created with
  [`brms::prior()`](https://paulbuerkner.com/brms/reference/set_prior.html))
  with names typically including "vague", "weakreg", and "informative"

- measure:

  Character string specifying the effect measure. One of "OR" (odds
  ratio), "RR" (risk ratio), "HR" (hazard ratio), "IRR" (incidence rate
  ratio), "MD" (mean difference), or "SMD" (standardized mean
  difference)

- study_var:

  Optional. Name of the study identifier variable (unquoted). Required
  if `incl_pet_peese = TRUE` or `incl_mixture = TRUE`

- rob_var:

  Optional. Name of the risk of bias variable (unquoted). Should contain
  values like "Low", "Some concerns", "High"

- exclude_high_rob:

  Logical. If TRUE, performs sensitivity analysis excluding studies with
  high risk of bias. Default is FALSE

- incl_pet_peese:

  Logical. If TRUE, includes PET-PEESE bias adjustment analysis. Default
  is FALSE

- pet_peese_direction:

  Character string. Direction of expected bias for PET-PEESE. Either
  "negative" or "positive". Default is "negative"

- pet_peese_threshold:

  Numeric. Threshold for choosing between PET and PEESE (based on
  proportion of studies showing effect). Default is 0.10

- incl_mixture:

  Logical. If TRUE, includes robust mixture model analysis. Default is
  FALSE

- incl_bma:

  Logical. If TRUE, includes Bayesian model averaging results. Default
  is FALSE

- model_bma:

  Named list of RoBMA model fits (from `RoBMA_brms()`) to include in
  model averaging. Names should be "vague", "weakreg", or "informative"

- prob_reference:

  Character string. What reference to use for probability calculations.
  Either "null" (compare to null_value) or "null_range" (compare to
  null_range boundaries). Default is "null"

- null_value:

  Numeric. The null value for the effect measure. If NULL, uses default
  based on measure (1 for ratio measures, 0 for difference measures)

- null_range:

  Numeric vector of length 2 or single numeric value. Defines the range
  of practical equivalence around the null value.

- add_null_range:

  Logical. If TRUE, displays the null range region on the plot. Default
  is FALSE

- color_null_range:

  Character string. Color for the null range region. Default is
  "#77bb41" (green)

- label_control:

  Character string. Label for the control/reference group. Default is
  "Control"

- label_intervention:

  Character string. Label for the intervention group. Default is
  "Intervention"

- title:

  Optional character string. Main title for the plot

- subtitle:

  Optional character string. Subtitle for the plot

- xlim:

  Optional numeric vector of length 2. X-axis limits for the plot

- x_breaks:

  Optional numeric vector. Custom x-axis break points

- color_palette:

  Optional named character vector of colors for each section. Names
  should match section labels (e.g., "All studies", "Excluding High
  RoB", etc.)

- show_density:

  Logical. If TRUE, includes a density plot below the ECDF. Default is
  TRUE

- prior_to_plot:

  Character string. Which prior specification to plot in the ECDF. One
  of "vague", "weakreg", or "informative". Default is "weakreg"

- font:

  Optional character string. Font family to use for the plot

## Value

A ggplot object (or patchwork object if `show_density = TRUE`)

## Details

This function creates an ECDF plot that displays the cumulative
probability distributions of pooled effects across different sensitivity
analyses. The plot includes:

- Left y-axis: Probability that mean effect \< x

- Right y-axis: Probability that mean effect \> x

- Optional density curves below the ECDF

- Optional region of practical equivalence (ROPE) shading

The `prob_reference` argument controls how probabilities are calculated:

- "null": Probabilities are calculated relative to the null value (e.g.,
  OR = 1)

- "null_range": Probabilities are calculated relative to the null range
  boundaries (left y-axis shows P(effect \< lower bound), right y-axis
  shows P(effect \> upper bound))

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic ECDF plot
sensitivity_ecdf_plot(
  model = my_brms_model,
  data = my_data,
  priors = list(
    vague = prior("normal(0, 10)", class = "Intercept"),
    weakreg = prior("normal(0, 1)", class = "Intercept"),
    informative = prior("normal(0.5, 0.5)", class = "Intercept")
  ),
  measure = "OR",
  prior_to_plot = "weakreg"
)

# Advanced ECDF with null range probability reference
sensitivity_ecdf_plot(
  model = my_model,
  data = my_data,
  priors = my_priors,
  measure = "OR",
  study_var = study_id,
  rob_var = rob,
  exclude_high_rob = TRUE,
  incl_pet_peese = TRUE,
  incl_mixture = TRUE,
  incl_bma = TRUE,
  model_bma = my_robma_models,
  prior_to_plot = "weakreg",
  prob_reference = "null_range",
  null_range = c(0.9, 1.1),
  add_null_range = TRUE
)
} # }
```
