#' Validate inputs for bayes_forest
#'
#' @noRd
validate_inputs <- function(
    model,
    data,
    measure,
    studyvar,
    yearvar,
    subgroup,
    subgroupvar,
    sort_studies_by
) {
  # Accept any supported model class (not just brmsfit)
  supported <- c("brmsfit", "meta_stan", "bayesmeta",
                  "stanreg", "stanfit", "CmdStanMCMC",
                  "bayesma", "brma", "rjags")
  if (!inherits(model, supported)) {
    rlang::abort(paste0(
      "`model` must be one of: ", paste(supported, collapse = ", "),
      ". Got: '", class(model)[1], "'"
    ))
  }
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.")
  }
  rlang::arg_match(measure, c("OR", "HR", "RR", "IRR", "MD", "SMD"))
  rlang::arg_match(sort_studies_by, c("author", "year", "effect"))

  if (rlang::quo_is_null(rlang::enquo(studyvar))) {
    rlang::abort("`studyvar` must be provided (bare column name).")
  }

  studyvar_name <- rlang::as_name(rlang::ensym(studyvar))
  if (!studyvar_name %in% names(data)) {
    rlang::abort(paste0("Can't find column `", studyvar_name, "` in `data`."))
  }

  if (!rlang::quo_is_null(rlang::enquo(yearvar))) {
    year_name <- rlang::as_name(rlang::ensym(yearvar))
    if (!year_name %in% names(data)) {
      rlang::abort(paste0("Can't find column `", year_name, "` in `data`."))
    }
  }

  # Validate subgroup variable if subgroup = TRUE
  if (isTRUE(subgroup)) {
    if (rlang::quo_is_null(rlang::enquo(subgroupvar))) {
      # Default to "Subgroup" column
      if (!"Subgroup" %in% names(data)) {
        rlang::abort("`subgroup = TRUE` but no `subgroupvar` specified and no 'Subgroup' column found in data.")
      }
    } else {
      # Check if specified subgroupvar exists
      subgroup_var_name <- rlang::as_name(rlang::ensym(subgroupvar))
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
#' @param model A fitted model object (any supported backend)
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
  supported <- c("brmsfit", "meta_stan", "bayesmeta",
                  "stanreg", "stanfit", "CmdStanMCMC")
  if (!inherits(model, supported)) {
    rlang::abort(paste0(
      "`model` must be one of: ", paste(supported, collapse = ", "),
      ". Got: '", class(model)[1], "'"
    ))
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
#' Checks that priors is a named list of prior specifications, where each
#' element contains \code{name}, \code{mu_prior} (a \code{bf_prior} object),
#' and optionally \code{tau_prior} (a \code{bf_prior} object).
#'
#' @param priors Object to validate.
#'
#' @return Invisible TRUE if valid, otherwise throws an error
#'
#' @keywords internal
#' @noRd
validate_priors <- function(priors) {

  if (is.null(priors)) {
    return(invisible(TRUE))
  }

  if (!is.list(priors) || is.null(names(priors))) {
    rlang::abort(
      "`priors` must be a named list. See `?bf_normal` for the expected format."
    )
  }

  for (nm in names(priors)) {
    spec <- priors[[nm]]

    if (!is.list(spec)) {
      rlang::abort(paste0(
        "Prior '", nm, "' must be a list with `name`, `mu_prior`, ",
        "and optionally `tau_prior`."
      ))
    }

    if (is.null(spec$name)) {
      rlang::abort(paste0("Prior '", nm, "' is missing a `name` element."))
    }

    if (!is.null(spec$mu_prior) && !inherits(spec$mu_prior, "bf_prior")) {
      rlang::abort(paste0(
        "Prior '", nm, "': `mu_prior` must be a bf_prior object ",
        "(e.g. bf_normal(0, 1)). Got: '", class(spec$mu_prior)[1], "'"
      ))
    }

    if (!is.null(spec$tau_prior) && !inherits(spec$tau_prior, "bf_prior")) {
      rlang::abort(paste0(
        "Prior '", nm, "': `tau_prior` must be a bf_prior object ",
        "(e.g. bf_half_cauchy(0, 0.5)). Got: '", class(spec$tau_prior)[1], "'"
      ))
    }
  }

  invisible(TRUE)
}
