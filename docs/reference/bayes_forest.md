# Create Bayesian Forest Plot for Meta-Analysis

Creates a Bayesian forest plot for meta-analysis using a fitted Bayesian
model. Supports multiple backends (brms, MetaStan, bayesmeta, rstanarm,
rstan, cmdstanr, bayesma, RoBMA, R2jags) and both common-effect and
random-effect models.

## Usage

``` r
bayes_forest(
  model,
  data,
  measure,
  studyvar = NULL,
  yearvar = NULL,
  n_ctrl = NULL,
  n_int = NULL,
  event_ctrl = NULL,
  event_int = NULL,
  mean_ctrl = NULL,
  mean_int = NULL,
  sd_ctrl = NULL,
  sd_int = NULL,
  time_ctrl = NULL,
  time_int = NULL,
  sort_studies_by = "author",
  model_type = c("random_effect", "common_effect"),
  re_min_k = NULL,
  subgroup = FALSE,
  subgroupvar = NULL,
  sort_subgroup_by = "alphabetical",
  label_outcome = "Outcome",
  label_control = "Control",
  label_intervention = "Intervention",
  title = NULL,
  subtitle = NULL,
  title_align = "left",
  xlim = NULL,
  x_breaks = NULL,
  x_label = NULL,
  effect_label = NULL,
  add_rope = FALSE,
  rope_value = NULL,
  rope_color = "grey50",
  plot_output = c("density", "boxplot"),
  shrinkage_output = c("density", "pointinterval"),
  null_value = NULL,
  color_palette = NULL,
  color_study_posterior_null_left = "deepskyblue",
  color_study_posterior_null_right = "violet",
  color_study_posterior = "dodgerblue",
  color_study_posterior_outline = "blue",
  color_overall_posterior = "blue",
  color_shrinkage_pointinterval = "purple",
  color_shrinkage_outline = "purple",
  color_shrinkage_fill = NULL,
  split_color_by_null = FALSE,
  color_favours_control = "firebrick",
  color_favours_intervention = "dodgerblue",
  add_arm_labels = TRUE,
  reverse_arm_labels = FALSE,
  add_pred = FALSE,
  add_pred_subgroup = FALSE,
  pred_output = c("density", "pointinterval", "boxplot"),
  color_pred_posterior = "orange",
  color_pred_outline = "darkorange2",
  color_pred_pointinterval = "orange",
  plot_width = 4,
  add_rob = FALSE,
  rob_tool = "rob2",
  add_rob_legend = FALSE,
  exclude_high_rob = FALSE,
  font = NULL
)
```

## Arguments

- model:

  A fitted Bayesian meta-analysis model object. Supported classes:
  `brmsfit` (brms), `meta_stan` (MetaStan), `bayesmeta`, `stanreg`
  (rstanarm), `stanfit` (rstan), `CmdStanMCMC` (cmdstanr), `bayesma`
  (bayesma), `brma` (RoBMA), `rjags` (R2jags/JAGS).

- data:

  A data frame containing the study data used for the meta-analysis.

- measure:

  Character string specifying the effect measure. Must be one of: "OR"
  (Odds Ratio), "HR" (Hazard Ratio), "RR" (Risk Ratio), "IRR" (Incidence
  Rate Ratio), "MD" (Mean Difference), or "SMD" (Standardized Mean
  Difference).

- studyvar:

  Column name containing study identifiers/authors.

- yearvar:

  Column name containing publication years. Default is NULL.

- n_ctrl:

  Column name containing control group sample sizes.

- n_int:

  Column name containing intervention group sample sizes.

- event_ctrl:

  Column name containing control group event counts.

- event_int:

  Column name containing intervention group event counts.

- mean_ctrl:

  Column name containing control group means.

- mean_int:

  Column name containing intervention group means.

- sd_ctrl:

  Column name containing control group standard deviations.

- sd_int:

  Column name containing intervention group standard deviations.

- time_ctrl:

  Column name containing control group time periods.

- time_int:

  Column name containing intervention group time periods.

- sort_studies_by:

  Character string specifying how to sort studies. Options: "author"
  (default), "year", or "effect".

- model_type:

  Pooling type used to draw the plot: "random_effect" or
  "common_effect". If not supplied, it is detected from `model`.

- re_min_k:

  Optional minimum number of studies for a random-effect model.
  Between-study heterogeneity (tau) is poorly identified when there are
  few studies. If `re_min_k` is set and there are fewer studies than
  this, a random-effect `model_type` is switched to "common_effect" with
  a warning. `bayesma` models are refit as common-effect (overall and
  per subgroup); for other backends only the plot changes, as the pooled
  estimate still comes from the fitted random-effect model, so refit the
  model as common-effect yourself. Default is NULL (no minimum). An
  alternative is to keep `re_min_k = NULL` and fit the model with an
  informative prior on tau.

- subgroup:

  Logical indicating whether to create subgroup analysis. Default is
  FALSE.

- subgroupvar:

  Column name in `data` to use for subgroup analysis.

- sort_subgroup_by:

  Character string or vector specifying subgroup ordering.

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

  Character string specifying title alignment.

- xlim:

  Numeric vector of length 2 specifying x-axis limits. Default is NULL.

- x_breaks:

  Numeric vector specifying custom x-axis break points. Default is NULL.

- x_label:

  Character string for the x-axis label. Default is NULL
  (measure-specific label).

- effect_label:

  Character string used for the effect columns in the right-hand table.
  Default is NULL (uses `measure`).

- add_rope:

  Logical indicating whether to add ROPE. Default is FALSE.

- rope_value:

  Numeric vector of length 2 specifying ROPE range.

- rope_color:

  Color for ROPE shading.

- plot_output:

  Character string specifying how study and pooled effects are drawn:
  "density" (default) or "boxplot". Boxplots are horizontal: the box
  spans the interquartile range, the line marks the median, and the
  whiskers span the 95% interval. With "boxplot", `shrinkage_output`
  must be "pointinterval" (it is set automatically if not supplied).

- shrinkage_output:

  Character string specifying shrinkage visualization. Options:
  "density" (default) or "pointinterval". Ignored for common-effect
  models.

- null_value:

  Numeric value specifying x-axis value of the null line. Default is
  NULL.

- color_palette:

  Character vector of colors for the plot. Default is NULL.

- color_study_posterior_null_left:

  Color for left side of null study posteriors.

- color_study_posterior_null_right:

  Color for right side of null study posteriors.

- color_study_posterior:

  Color for study posterior densities.

- color_study_posterior_outline:

  Color for study posterior outlines.

- color_overall_posterior:

  Color for overall posterior.

- color_shrinkage_pointinterval:

  Color for point intervals.

- color_shrinkage_outline:

  Color for shrinkage plot outlines.

- color_shrinkage_fill:

  Color for shrinkage plot fill.

- split_color_by_null:

  Logical. Split posterior densities by null value. For boxplots, each
  box is coloured by which side of the null its median falls.

- color_favours_control:

  Colour for density regions favouring the control group.

- color_favours_intervention:

  Colour for density regions favouring the intervention group.

- add_arm_labels:

  Logical. Whether to show "Favours Control/Intervention" labels.
  Default is TRUE.

- reverse_arm_labels:

  Logical. Whether to swap the arm label positions. Default is FALSE.

- add_pred:

  Logical. Whether to add a predictive distribution for a new study.
  Default is FALSE. Only applicable to random-effect models.

- add_pred_subgroup:

  Logical. Whether to add prediction to each subgroup. Default is FALSE.

- pred_output:

  Character string specifying prediction visualization. Options:
  "density" (default), "pointinterval", or "boxplot". When
  `plot_output = "boxplot"`, `pred_output` defaults to "boxplot" too
  unless explicitly set.

- color_pred_posterior:

  Color for prediction density fill. Default is "orange".

- color_pred_outline:

  Color for prediction density outline. Default is "darkorange2".

- color_pred_pointinterval:

  Color for prediction point interval. Default is "orange".

- plot_width:

  Numeric value specifying the relative width of the plot component.
  Default is 4.

- add_rob:

  Logical indicating whether to add Risk of Bias assessment. Default is
  FALSE.

- rob_tool:

  Risk of bias tool: "rob2" (default), "robins_i", "robins_e",
  "quadas2", "clarity", "rob_nma", "probast", or "probast_ai". See
  [`rob_plot()`](https://blmoran.github.io/bayesfoRest/reference/rob_plot.md)
  for the domain columns each tool expects.

- add_rob_legend:

  Logical indicating whether to add RoB legend. Default is FALSE.

- exclude_high_rob:

  Logical indicating whether to exclude studies with a high overall risk
  of bias ("High", "Serious", "Critical", "Very high" or "Definitely
  high") and refit the model. Requires an overall RoB column.

- font:

  Character string specifying the font family. Default is NULL.

## Value

A patchwork object containing the complete forest plot.
