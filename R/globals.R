# Global variable declarations to silence R CMD check NOTEs
# These are column names used in dplyr/tidybayes pipelines via NSE

utils::globalVariables(c(
  ".draw", ".grp", ".last_in_group", ".lower", ".lower_tau", ".upper", ".upper_tau",
  ":=",
  "Author", "Author_label", "Author_ordered", "Author_pooled",
  "Check", "Code",
  "D1", "Description",
  "N_ctrl", "N_int",
  "Overall",
  "Subgroup",
  "Year",
  "band",
  "category", "conflict", "contraction", "ctrl_mean_sd",
  "effect_draws", "est", "estimate",
  "flag",
  "hi",
  "int_mean_sd", "interpretation",
  "l95", "level", "level_num", "lo",
  "max_level", "mean_effect", "median", "mu", "mu_hi", "mu_lo", "mu_median",
  "mu_prior_unicode",
  "panel", "parameter", "plot_row", "pr_benefit", "pred", "prior",
  "prior_label", "prior_percentile", "prior_sd",
  "ratio", "re",
  "section_label", "sei", "side", "spread_df", "study_count",
  "subgroup.forest.summary", "subgroup_model",
  "tau", "tau_hi", "tau_lo", "tau_median", "tau_prior", "tau_prior_pctile",
  "tau_prior_sd", "tau_prior_unicode", "theta", "theta_study",
  "u95", "unweighted_effect",
  "val", "variable", "vi",
  "weighted_effect",
  "x", "x_studies", "x_val", "xdist",
  "y", "yi", "yint",
  ".above_null", ".fill", ".overall", "n", "se",
  "xlower", "xmax", "xmiddle", "xmin", "xupper"
))

#' @importFrom rlang .data
#' @rawNamespace if (getRversion() < "4.3.0") importFrom("S7", "@")
NULL

.onLoad <- function(libname, pkgname) {
  S7::methods_register()
}

#' @importFrom stats dexp dnorm median update
#' @importFrom utils capture.output head
NULL
