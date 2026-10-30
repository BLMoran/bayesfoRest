# Global variable declarations to silence R CMD check NOTEs
# These are column names used in dplyr/tidybayes pipelines via NSE

utils::globalVariables(c(
  ".draw", ".last_in_group", ".lower", ".lower_tau", ".upper", ".upper_tau",
  ":=",
  "Author", "Author_ordered", "Author_pooled",
  "Code",
  "D1", "Description",
  "Overall",
  "Subgroup",
  "Year",
  "category", "ctrl_mean_sd",
  "effect_draws", "estimate",
  "int_mean_sd",
  "l95", "level_num",
  "max_level", "mean_effect", "median", "mu", "mu_prior_unicode",
  "plot_row", "pred", "prior", "prior_label",
  "section_label", "spread_df", "study_count",
  "subgroup.forest.summary", "subgroup_model",
  "tau", "tau_prior_unicode", "theta_study",
  "u95", "unweighted_effect",
  "vi",
  "weighted_effect",
  "x", "x_studies", "xdist",
  "yi",
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
