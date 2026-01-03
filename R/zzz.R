
# Global variables used in tidyverse pipelines
utils::globalVariables(c(
  "Author", "Author_ordered", "Author_pooled", "Year", 
  "Subgroup", "Overall", "D1", "D2", "D3", "D4", "D5",
  "yi", "vi", "b_Intercept", "sd_Author__Intercept", "r_Author",
  ".lower", ".upper", "weighted_effect", "unweighted_effect",
  "int_mean_sd", "ctrl_mean_sd", "effect_draws", "study_count",
  "subgroup_model", "subgroup.forest.summary", "spread_df",
  "mean_effect", "null_line_to_use", "x_studies", "xdist",
  "Code", "Description", "Event_Control", "Event_Intervention",
  "N_Control", "N_Intervention", "N_Total",
  "Outcome_Control_No", "Outcome_Intervention_No", "x",
  # overall_plot
  "mu", "tau", "category", "density", "coef",
  
  # sensitivity functions
  "section_label", "prior", "prior_label", "plot_row", "level_num", "max_level",
  "mu_prior_unicode", "tau_prior_unicode", ". last_in_group",
  "median", "l95", "u95", "estimate", "post_prob",
  
  # other
  "normal", "formula", "update", "as. formula"
))

