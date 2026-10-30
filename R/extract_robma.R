#' Extract Posterior Draws from RoBMA for Sensitivity Analysis
#'
#' @description
#' Extracts model-averaged posterior draws for the overall effect (mu) from a
#' fitted RoBMA object, optionally filtering by model component type. This
#' allows the sensitivity plot to display separate rows for PET-PEESE,
#' selection (weight) models, and full RoBMA (conditional vs model-averaged).
#'
#' @param robma_fit A fitted RoBMA object (from \code{RoBMA::RoBMA()}).
#' @param type Character string specifying which component draws to extract.
#'   One of:
#'   \describe{
#'     \item{"model_averaged"}{Model-averaged mu across all models, including
#'       null effect models (spike at zero). Default.}
#'     \item{"conditional"}{Mu conditional on effect being present (non-null
#'       effect models only).}
#'     \item{"pet_peese"}{Mu conditional on PET and/or PEESE publication bias
#'       component being active.}
#'     \item{"selection"}{Mu conditional on selection (weight function) component
#'       being active.}
#'   }
#' @param measure Character string specifying the effect measure (e.g. "OR",
#'   "RR", "SMD"). Used to determine whether to exponentiate draws.
#'
#' @return A numeric vector of posterior draws on the outcome scale (i.e.
#'   exponentiated for ratio measures).
#'
#' @details
#' For \code{"pet_peese"} and \code{"selection"} types, the function uses
#' \code{BayesTools::as_mixed_posteriors()} to extract draws from the
#' product-space model conditional on the relevant bias component indicator.
#' If the requested component is not present in the model, the function
#' returns \code{NULL} with a warning.
#'
#' @examples
#' \dontrun{
#' # Extract model-averaged draws
#' draws_ma <- extract_robma_sensitivity_draws(fit, "model_averaged", "OR")
#'
#' # Extract PET-PEESE bias-corrected draws
#' draws_pp <- extract_robma_sensitivity_draws(fit, "pet_peese", "OR")
#' }
#'
#' @keywords internal
#' @noRd
extract_robma_sensitivity_draws <- function(
    robma_fit,
    type = c("model_averaged", "conditional", "pet_peese", "selection"),
    measure
) {

  type <- match.arg(type)

  if (!requireNamespace("RoBMA", quietly = TRUE)) {
    rlang::abort("Package 'RoBMA' is required for RoBMA models.")
  }

  draws <- switch(
    type,

    model_averaged = {
      tryCatch(
        as.numeric(RoBMA::pooled_effect(robma_fit, conditional = FALSE)),
        error = function(e) {
          rlang::warn(paste0("RoBMA model-averaged extraction failed: ", e$message))
          NULL
        }
      )
    },

    conditional = {
      tryCatch(
        as.numeric(RoBMA::pooled_effect(robma_fit, conditional = TRUE)),
        error = function(e) {
          rlang::warn(paste0("RoBMA conditional extraction failed: ", e$message))
          NULL
        }
      )
    },

    pet_peese = {
      .extract_robma_component_draws(robma_fit, component = "PETPEESE")
    },

    selection = {
      .extract_robma_component_draws(robma_fit, component = "omega")
    }
  )

  if (is.null(draws)) return(NULL)

  if (measure %in% c("OR", "RR", "HR", "IRR")) {
    draws <- exp(draws)
  }

  draws
}


#' Extract Component-Conditional Draws from RoBMA Product-Space Model
#'
#' @description
#' Internal helper that extracts mu posterior draws from the product-space JAGS
#' fit, conditional on a specific publication-bias component being active.
#' Uses the indicator variable approach in \code{BayesTools::as_mixed_posteriors()}.
#'
#' @param robma_fit A fitted RoBMA object.
#' @param component Character. Conditioning component passed to
#'   \code{BayesTools::as_mixed_posteriors()} via its \code{conditional}
#'   argument. One of \code{"PETPEESE"}, \code{"PET"}, \code{"PEESE"},
#'   \code{"omega"}.
#'
#' @return Numeric vector of posterior draws, or NULL if the component
#'   is not present.
#'
#' @keywords internal
#' @noRd
.extract_robma_component_draws <- function(robma_fit, component) {

  if (!requireNamespace("BayesTools", quietly = TRUE)) {
    rlang::abort("Package 'BayesTools' is required for component-level RoBMA extraction.")
  }

  fit <- robma_fit[["fit"]]
  if (is.null(fit) || length(fit) == 0) {
    rlang::warn("RoBMA object does not contain a valid fit.")
    return(NULL)
  }

  mu_par <- .robma_mu_parameter(robma_fit)

  tryCatch({
    mixed <- BayesTools::as_mixed_posteriors(
      model      = fit,
      parameters = mu_par,
      conditional = component
    )
    as.numeric(mixed[[mu_par]])
  }, error = function(e) {
    rlang::warn(paste0(
      "Failed to extract ", component, "-conditional draws: ", e$message
    ))
    NULL
  })
}


#' Detect the mu Parameter Name in a RoBMA Fit
#'
#' @param robma_fit A fitted RoBMA object.
#' @return Character scalar: \code{"mu"} or \code{"mu_intercept"}.
#' @keywords internal
#' @noRd
.robma_mu_parameter <- function(robma_fit) {
  prior_list <- attr(robma_fit[["fit"]], "prior_list")
  if ("mu" %in% names(prior_list)) "mu" else "mu_intercept"
}


#' Convert RoBMA Draws to Sensitivity Format
#'
#' @description
#' Wraps \code{extract_robma_sensitivity_draws()} and returns draws in the
#' tibble format expected by the sensitivity plot pipeline.
#'
#' @param robma_fit A fitted RoBMA object
#' @param type Character. Model component type (see
#'   \code{extract_robma_sensitivity_draws}).
#' @param measure Character. Effect measure.
#' @param prior Character. Prior name ("vague", "weakreg", "informative").
#' @param prior_label Character. Display label for the prior.
#' @param section_label Character. Section label for the sensitivity plot row group.
#'
#' @return A tibble with columns: x, prior, prior_label, section_label
#'
#' @keywords internal
#' @noRd
robma_to_sensitivity_draws <- function(
    robma_fit,
    type    = "model_averaged",
    measure,
    prior,
    prior_label,
    section_label
) {

  draws <- extract_robma_sensitivity_draws(
    robma_fit = robma_fit,
    type      = type,
    measure   = measure
  )

  if (is.null(draws)) {
    return(tibble::tibble(
      x             = numeric(),
      prior         = character(),
      prior_label   = character(),
      section_label = character()
    ))
  }

  tibble::tibble(
    x             = draws,
    prior         = prior,
    prior_label   = prior_label,
    section_label = section_label
  )
}


#' Build All RoBMA Sensitivity Draws
#'
#' @description
#' Given a named list of RoBMA fits (one per prior level), extracts draws for
#' each requested component type and returns a combined tibble ready for the
#' sensitivity plot.
#'
#' @param robma_fits Named list of fitted RoBMA objects. Names should be prior
#'   keys (e.g. "vague", "weakreg", "informative").
#' @param measure Character. Effect measure.
#' @param types Character vector of component types to extract. Defaults to
#'   \code{c("pet_peese", "selection", "conditional", "model_averaged")}.
#' @param section_labels Named character vector mapping type names to display
#'   labels. Defaults provided.
#'
#' @return A tibble with columns: x, prior, prior_label, section_label
#'
#' @keywords internal
#' @noRd
build_robma_sensitivity_draws <- function(
    robma_fits,
    measure,
    types = c("pet_peese", "selection", "conditional", "model_averaged"),
    section_labels = c(
      pet_peese      = "PET-PEESE",
      selection      = "Selection Model (Weight)",
      conditional    = "RoBMA (Conditional)",
      model_averaged = "RoBMA (Model-Averaged)"
    )
) {

  prior_labels <- c(
    vague       = "Vague",
    weakreg     = "Weakly Regularising",
    informative = "Informative",
    single      = ""
  )

  purrr::map(types, function(type) {
    purrr::imap(robma_fits, function(fit, prior_name) {
      robma_to_sensitivity_draws(
        robma_fit     = fit,
        type          = type,
        measure       = measure,
        prior         = prior_name,
        prior_label   = prior_labels[[prior_name]] %||% prior_name,
        section_label = section_labels[[type]] %||% type
      )
    }) |> purrr::list_rbind()
  }) |> purrr::list_rbind()
}
