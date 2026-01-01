#' Update Model with New Prior
#'
#' @description
#' Updates a brms model with new prior specifications while keeping the same
#' data and model structure.
#'
#' @param model Original brmsfit object
#' @param data Data frame for the model
#' @param new_prior New prior specification(s)
#'
#' @return Updated brmsfit object
#'
#' @keywords internal
#' @noRd
update_model_prior <- function(model, data, new_prior) {
  update(
    model,
    newdata = data,
    prior = new_prior,
    refresh = 0,
    recompile = FALSE
  )
}

#' Perform Prior Sensitivity Analysis
#'
#' @description
#' Re-fits a model with different prior specifications and extracts posterior draws
#' for sensitivity analysis.
#'
#' @param model Original brmsfit object
#' @param data Data frame containing the study data
#' @param priors Named list of prior specifications
#' @param measure Effect measure type
#' @param section_label Label for this section of the analysis
#' @param prior_labels Named vector mapping prior names to display labels
#'
#' @return A tibble with columns: section_label, prior, prior_label, and x (draws)
#'
#' @keywords internal
#' @noRd
update_prior_sensitivity <- function(model,
                                     data,
                                     priors,
                                     measure,
                                     section_label,
                                     prior_labels) {
  
  purrr::imap_dfr(priors, function(prior_obj, prior_name) {
    
    updated_model <- update_model_prior(
      model = model,
      data = data,
      new_prior = prior_obj
    )
    
    tibble::tibble(
      section_label = section_label,
      prior = prior_name,
      prior_label = prior_labels[[prior_name]],
      x = extract_intercept_draws(updated_model, measure)
    )
  })
}

#' Summarize Sensitivity Analysis Posteriors
#'
#' @description
#' Computes summary statistics and probabilities from posterior draws across
#' all sensitivity analyses.
#'
#' @param draws A tibble containing posterior draws from sensitivity analyses
#' @param null_value The null value for hypothesis testing
#' @param null_range Optional range for region of practical equivalence
#'
#' @return A tibble with summary statistics including median, 95% CrI, and probabilities
#'
#' @keywords internal
#' @noRd
summarise_sensitivity_posteriors <- function(
    draws,
    null_value,
    null_range
) {
  
  section_levels <- unique(draws$section_label)
  
  prior_levels <- c("vague", "weakreg", "informative")
  
  draws |>
    dplyr::mutate(
      section_label = factor(section_label, levels = section_levels),
      prior         = factor(prior, levels = prior_levels)
    ) |>
    dplyr::group_by(section_label, prior, prior_label) |>
    dplyr::summarise(
      median = stats::median(x),
      l95    = stats::quantile(x, 0.025),
      u95    = stats::quantile(x, 0.975),
      
      pr_benefit = round(mean(x < null_value) * 100, 1),
      pr_no_benefit = round(mean(x > null_value) * 100, 1),
      
      pr_benefit_null_range =
        round(mean(x < null_range[1]) * 100, 1),
      pr_no_benefit_null_range =
        round(mean(x > null_range[2]) * 100, 1),
      
      .groups = "drop"
    ) |>
    dplyr::arrange(section_label, prior)
}

#' Convert RoBMA Results to Sensitivity Draws
#'
#' @description
#' Converts model-averaged posterior draws from RoBMA to the format
#' needed for sensitivity analysis plotting.
#'
#' @param robma_fit A RoBMA_brms fit object
#' @param prior Prior name ("vague", "weakreg", or "informative")
#' @param prior_label Display label for the prior
#' @param section_label Label for this analysis section
#'
#' @return A tibble with columns: x, prior, prior_label, section_label
#'
#' @keywords internal
#' @noRd
robma_to_sensitivity_draws <- function(robma_fit,
                                       prior,
                                       prior_label,
                                       section_label = "Model-averaged") {
  x <- robma_fit$ma_posterior
  
  if (robma_fit$measure %in% c("OR", "RR", "HR", "IRR")) {
    x <- exp(x)
  }
  
  tibble::tibble(
    x = x,
    prior = prior,
    prior_label = prior_label,
    section_label = section_label
  )
}
