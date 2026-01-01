#' Validate inputs for bayes_forest
#'
#' @noRd
validate_inputs <- function(
    model,
    data,
    measure,
    studyvar,
    year,
    subgroup,
    subgroup_var,
    sort_studies_by,
    shrinkage_output
) {
  if (!inherits(model, "brmsfit")) {
    rlang::abort("`model` must be a brmsfit object.")
  }
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.")
  }
  rlang::arg_match(measure, c("OR", "HR", "RR", "IRR", "MD", "SMD"))
  rlang::arg_match(sort_studies_by, c("author", "year", "effect"))
  rlang::arg_match(shrinkage_output, c("density", "pointinterval"))
  
  if (rlang::quo_is_null(rlang::enquo(studyvar))) {
    rlang::abort("`studyvar` must be provided (bare column name).")
  }
  
  studyvar_name <- rlang::as_name(rlang::ensym(studyvar))
  if (!studyvar_name %in% names(data)) {
    rlang::abort(paste0("Can't find column `", studyvar_name, "` in `data`."))
  }
  
  if (!rlang::quo_is_null(rlang::enquo(year))) {
    year_name <- rlang::as_name(rlang::ensym(year))
    if (!year_name %in% names(data)) {
      rlang::abort(paste0("Can't find column `", year_name, "` in `data`."))
    }
  }
  
  # Validate subgroup variable if subgroup = TRUE
  if (isTRUE(subgroup)) {
    if (rlang::quo_is_null(rlang::enquo(subgroup_var))) {
      # Default to "Subgroup" column
      if (!"Subgroup" %in% names(data)) {
        rlang::abort("`subgroup = TRUE` but no `subgroup_var` specified and no 'Subgroup' column found in data.")
      }
    } else {
      # Check if specified subgroup_var exists
      subgroup_var_name <- rlang::as_name(rlang::ensym(subgroup_var))
      if (!subgroup_var_name %in% names(data)) {
        rlang::abort(paste0("Can't find subgroup column `", subgroup_var_name, "` in `data`."))
      }
    }
  }
  
  invisible(TRUE)
}

#' Validate Inputs for Sensitivity Plot
#'
#' @description
#' Internal validation function that checks the validity of core inputs to the
#' sensitivity plot function.
#'
#' @param model A brmsfit object to validate
#' @param data A data frame to validate
#' @param measure Character string specifying the effect measure
#'
#' @return Invisible TRUE if validation passes, otherwise throws an error
#'
#' @keywords internal
#' @noRd
validate_inputs_sens_plot <- function(
    model,
    data,
    priors,
    measure
) {
  
  if (!inherits(model, "brmsfit")) {
    rlang::abort("`model` must be a brmsfit object.")
  }
  
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.")
  }
  
  rlang::arg_match(measure, c("OR", "HR", "RR", "IRR", "MD", "SMD"))
  
  validate_priors(priors)
  
  invisible(TRUE)
}

#' Validate Priors Object
#'
#' @description
#' Checks that priors is either a single brms prior or a list of brms priors.
#'
#' @param priors Object to validate as brms prior(s)
#'
#' @return Invisible TRUE if valid, otherwise throws an error
#'
#' @keywords internal
#' @noRd
validate_priors <- function(priors) {
  
  if (is.null(priors)) {
    return(invisible(TRUE))
  }
  
  if (inherits(priors, "brmsprior")) {
    return(invisible(TRUE))
  }
  
  if (is.list(priors) && all(vapply(priors, inherits, logical(1), "brmsprior"))) {
    return(invisible(TRUE))
  }
  
  rlang::abort(
    "`priors` must be a brms prior or a list of brms priors created with brms::prior()."
  )
}