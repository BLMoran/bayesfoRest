# Generate Sensitivity Analysis Plot for Bayesian Meta-Analysis

Creates a comprehensive sensitivity analysis visualization showing how
meta-analytic estimates vary across different prior specifications and
analysis approaches. The plot combines forest-plot style density
distributions with tabulated results.

## Usage

``` r
sensitivity_plot(
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
  add_probs = FALSE,
  null_value = NULL,
  null_range = NULL,
  add_null_range = FALSE,
  color_null_range = "#77bb41",
  label_control = "Control",
  label_intervention = "Intervention",
  title = NULL,
  subtitle = NULL,
  title_align = "left",
  xlim = NULL,
  x_breaks = NULL,
  color_palette = NULL,
  color_overall_posterior = "dodgerblue",
  color_overall_posterior_outline = "blue",
  split_color_by_null = FALSE,
  color_favours_control = "firebrick",
  color_favours_intervention = "dodgerblue",
  plot_width = 4,
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

- add_probs:

  Logical. If TRUE, adds probability columns to the results table.
  Default is FALSE

- null_value:

  Numeric. The null value for the effect measure. If NULL, uses default
  based on measure (1 for ratio measures, 0 for difference measures)

- null_range:

  Numeric vector of length 2 or single numeric value. Defines the range
  of practical equivalence around the null value. If single value,
  creates symmetric range around null_value

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

- title_align:

  Character string. Alignment for title and subtitle. One of "left",
  "center", or "right". Default is "left"

- xlim:

  Optional numeric vector of length 2. X-axis limits for the density
  plot

- x_breaks:

  Optional numeric vector. Custom x-axis break points

- color_palette:

  Optional color palette for the plot

- color_overall_posterior:

  Character string. Fill color for posterior densities. Default is
  "dodgerblue"

- color_overall_posterior_outline:

  Character string. Outline color for posterior densities. Default is
  "blue"

- split_color_by_null:

  Logical. If TRUE, colors densities based on which side of null they
  favor. Default is FALSE

- color_favours_control:

  Character string. Color when effect favors control. Default is
  "firebrick"

- color_favours_intervention:

  Character string. Color when effect favors intervention. Default is
  "dodgerblue"

- plot_width:

  Numeric. Width ratio for the density plot section. Default is 4

- font:

  Optional character string. Font family to use for the plot

## Value

A patchwork object combining the sensitivity analysis tables and density
plots

## Details

This function performs comprehensive sensitivity analyses including:

- Prior sensitivity: Re-fits the model with different prior
  specifications

- Risk of bias sensitivity: Excludes high risk of bias studies

- PET-PEESE: Publication bias adjustment using precision-based methods

- Mixture models: Robust analysis using mixture distributions

- Bayesian model averaging: Combines results across model specifications

The output combines three elements using patchwork:

- Left table: Prior specifications

- Center plot: Density distributions of posterior estimates

- Right table: Numerical summaries and probabilities

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic sensitivity plot
sensitivity_plot(
  model = my_brms_model,
  data = my_data,
  priors = list(
    vague = prior("normal(0, 10)", class = "Intercept"),
    weakreg = prior("normal(0, 1)", class = "Intercept"),
    informative = prior("normal(0.5, 0.5)", class = "Intercept")
  ),
  measure = "SMD"
)

# Advanced sensitivity with all options
sensitivity_plot(
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
  add_probs = TRUE,
  null_range = c(0.9, 1.1),
  add_null_range = TRUE
)
} # }
```
