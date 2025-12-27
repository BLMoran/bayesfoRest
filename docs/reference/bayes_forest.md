# Create Bayesian Forest Plot for Meta-Analysis

This function creates a Bayesian forest plot for meta-analysis using
data prepared with metafor and a fitted brms model. It supports various
effect measures and can handle both single and subgroup analyses with
customizable visualization options.

## Usage

``` r
bayes_forest(
  model,
  data,
  measure,
  studyvar = NULL,
  year = NULL,
  c_n = NULL,
  i_n = NULL,
  c_event = NULL,
  i_event = NULL,
  c_mean = NULL,
  i_mean = NULL,
  c_sd = NULL,
  i_sd = NULL,
  c_time = NULL,
  i_time = NULL,
  sort_studies_by = "author",
  subgroup = FALSE,
  subgroup_var = NULL,
  sort_subgroup_by = "alphabetical",
  label_outcome = "Outcome",
  label_control = "Control",
  label_intervention = "Intervention",
  title = NULL,
  subtitle = NULL,
  title_align = "left",
  xlim = NULL,
  x_breaks = NULL,
  add_rope = FALSE,
  rope_value = NULL,
  rope_color = "grey50",
  shrinkage_output = "density",
  null_value = NULL,
  color_palette = NULL,
  color_study_posterior_null_left = "deepskyblue",
  color_study_posterior_null_right = "violet",
  color_study_posterior = "dodgerblue",
  color_study_posterior_outline = "blue",
  color_overall_posterior = "blue",
  color_pointinterval = "purple",
  color_shrinkage_outline = "purple",
  color_shrinkage_fill = NULL,
  split_color_by_null = FALSE,
  color_favours_control = "firebrick",
  color_favours_intervention = "dodgerblue",
  plot_width = 4,
  add_rob = FALSE,
  rob_tool = c("rob2", "robins_i", "quadas2", "robins_e"),
  add_rob_legend = FALSE,
  exclude_high_rob = FALSE,
  font = NULL
)
```

## Arguments

- model:

  A fitted brms model object (class 'brmsfit') containing the Bayesian
  meta-analysis results.

- data:

  A data frame containing the study data used for the meta-analysis.

- measure:

  Character string specifying the effect measure. Must be one of: "OR"
  (Odds Ratio), "HR" (Hazard Ratio), "RR" (Risk Ratio), "IRR" (Incidence
  Rate Ratio), "MD" (Mean Difference), or "SMD" (Standardized Mean
  Difference).

- studyvar:

  Column name containing study identifiers/authors. Default is NULL.

- year:

  Column name containing publication years. Default is NULL.

- c_n:

  Column name containing control group sample sizes. Required for OR,
  RR, MD, SMD.

- i_n:

  Column name containing intervention group sample sizes. Required for
  OR, RR, MD, SMD.

- c_event:

  Column name containing control group event counts. Required for OR,
  RR, IRR.

- i_event:

  Column name containing intervention group event counts. Required for
  OR, RR, IRR.

- c_mean:

  Column name containing control group means. Required for MD, SMD.

- i_mean:

  Column name containing intervention group means. Required for MD, SMD.

- c_sd:

  Column name containing control group standard deviations. Required for
  MD, SMD.

- i_sd:

  Column name containing intervention group standard deviations.
  Required for MD, SMD.

- c_time:

  Column name containing control group time periods. Required for IRR.

- i_time:

  Column name containing intervention group time periods. Required for
  IRR.

- sort_studies_by:

  Character string specifying how to sort studies. Options: "author"
  (default), "year", or "effect".

- subgroup:

  Logical indicating whether to create subgroup analysis. Default is
  FALSE.

- subgroup_var:

  Character string. Name of the variable in data to use for subgroup
  analysis.

- sort_subgroup_by:

  Character string or vector specifying subgroup ordering. Options:
  "alphabetical" (default), "effect", or custom character vector of
  subgroup names.

- label_outcome:

  Character string for outcome label. Default is "Outcome".

- label_control:

  Character string for control group label. Default is "Control".

- label_intervention:

  Character string for intervention group label. Default is
  "Intervention".

- title:

  Character string for the plot title. Default is NULL (no title).

- subtitle:

  Character string for the plot subtitle. Default is NULL (no subtitle).

- title_align:

  Character string specifying title alignment. Options: "left"
  (default), "center"/"centre", "right".

- xlim:

  Numeric vector of length 2 specifying x-axis limits. Default is NULL
  (auto-scaled).

- x_breaks:

  Numeric vector specifying custom x-axis break points. Default is NULL
  (uses measure-specific defaults).

- add_rope:

  Logical indicating whether to add ROPE (Region of Practical
  Equivalence). Default is FALSE.

- rope_value:

  Numeric vector of length 2 specifying ROPE range, or single value for
  symmetric range around null. Default is NULL (uses Kruschke's
  recommendations: OR/HR/RR/IRR = c(0.9, 1.1), SMD = c(-0.1, 0.1), MD
  requires specification).

- rope_color:

  Color for ROPE shading. Default is transparent grey.

- shrinkage_output:

  Character string specifying shrinkage visualization. Options:
  "density" (default) or "pointinterval".

- null_value:

  Numeric value specifying x-axis value of the null line. Default is
  NULL (uses measure-specific defaults).

- color_palette:

  Character vector of colors for the plot. Default is NULL.

- color_study_posterior_null_left:

  Color for left side of null study posteriors. Default is
  "deepskyblue".

- color_study_posterior_null_right:

  Color for right side of null study posteriors. Default is "violet".

- color_study_posterior:

  Color for study posterior densities. Default is "dodgerblue".

- color_study_posterior_outline:

  Color for study posterior outlines. Default is "blue".

- color_overall_posterior:

  Color for overall posterior. Default is "blue".

- color_pointinterval:

  Color for point intervals. Default is "purple".

- color_shrinkage_outline:

  Color for shrinkage plot outlines. Default is "purple".

- color_shrinkage_fill:

  Color for shrinkage plot fill. Default is NULL.

- split_color_by_null:

  Logical. If TRUE, posterior densities are split and coloured based on
  whether values fall above or below the null value.

- color_favours_control:

  Colour used for density regions favouring the control group when
  `split_color_by_null = TRUE`.

- color_favours_intervention:

  Colour used for density regions favouring the intervention group when
  `split_color_by_null = TRUE`.

- plot_width:

  Numeric value specifying the relative width of the plot component.
  Default is 4.

- add_rob:

  Logical indicating whether to add Risk of Bias assessment. Default is
  FALSE.

- rob_tool:

  Character string specifying RoB tool. Options: "rob2" (default).

- add_rob_legend:

  Logical indicating whether to add RoB legend. Default is FALSE.

- exclude_high_rob:

  Logical indicating whether to exclude high risk of bias studies and
  refit the model. Default is FALSE.

- font:

  Character string specifying the font family to use throughout the
  plot. Default is NULL (uses system defaults). Common options include
  "Arial", "Times New Roman", "Helvetica", "Georgia", etc. The font must
  be available on your system.

## Value

A patchwork object containing the complete forest plot with study
information table, density plots, and effect size table.

## Details

The function performs several key operations:

- Validates input parameters and data structure

- Renames data columns based on the specified measure type

- Updates the brms model if necessary

- Handles Risk of Bias data if specified

- Creates either single or subgroup forest plots

- Combines tables and plots using patchwork

For Risk of Bias functionality, the data must contain columns named
"D1", "D2", "D3", and "Overall".

When `exclude_high_rob = TRUE`, studies with "High" overall risk of bias
are removed and the model is refitted on the filtered data.

## See also

[`brm`](https://paulbuerkner.com/brms/reference/brm.html) for fitting
Bayesian models,
[`metafor`](https://wviechtb.github.io/metafor/reference/metafor-package.html)
for creating meta analysis models,
[`gt`](https://gt.rstudio.com/reference/gt.html) for making great
tables,
[`patchwork`](https://patchwork.data-imaginist.com/reference/patchwork-package.html)
for combining plots

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with odds ratio
forest_plot <- bayes_forest(
  model = my_brms_model,
  data = my_data,
  measure = "OR",
  studyvar = author,
  year = year,
  c_n = control_n,
  i_n = intervention_n,
  c_event = control_events,
  i_event = intervention_events
)

# Subgroup analysis with mean difference
forest_plot_subgroup <- bayes_forest(
  model = my_brms_model,
  data = my_data,
  measure = "MD",
  studyvar = author,
  c_n = control_n,
  i_n = intervention_n,
  c_mean = control_mean,
  i_mean = intervention_mean,
  c_sd = control_sd,
  i_sd = intervention_sd,
  subgroup = TRUE,
  sort_subgroup_by = "effect"
)
} # }
```
