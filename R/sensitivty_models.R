#' Run Meta-Analysis Sensitivity
#'
#' @description
#' Dispatcher function that runs the appropriate sensitivity analysis
#' based on the specification provided.
#'
#' @param spec List containing analysis specifications
#' @param priors Named list of prior specifications
#' @param measure Effect measure type
#' @param prior_labels Named vector of prior labels
#'
#' @return A tibble with posterior draws and labels
#'
#' @keywords internal
#' @noRd
run_ma_sensitivity <- function(spec, priors, measure, prior_labels) {
  
  if (spec$id == "pet_peese" && !is.null(spec$model)) {
    
    # PET-PEESE handling
    pet_peese_to_sensitivity_draws(
      model         = spec$model,
      data          = spec$data,
      priors        = priors,
      measure       = measure,
      section_label = spec$section_label,
      prior_labels  = prior_labels,
      study_var     = spec$study_var,
      direction     = spec$pet_peese_direction %||% "negative",
      threshold     = spec$pet_peese_threshold %||% 0.10
    )
    
  } else if (spec$id == "mixture" && !is.null(spec$model)) {
    
    # Robust mixture sensitivity
    mixture_model_to_sensitivity_draws(
      model         = spec$model,
      data          = spec$data,
      priors        = priors,
      measure       = measure,
      section_label = spec$section_label,
      prior_labels  = prior_labels,
      study_var     = spec$study_var
    )
    
  } else {
    
    # Standard prior sensitivity
    update_prior_sensitivity(
      model         = spec$model,
      data          = spec$data,
      priors        = priors,
      measure       = measure,
      section_label = spec$section_label,
      prior_labels  = prior_labels
    )
  }
}

#' Mixture Model Sensitivity Analysis
#'
#' @description
#' Performs sensitivity analysis using robust mixture models that combine
#' a standard component with a robust component to handle outliers.
#'
#' @param model Original brmsfit object
#' @param data Data frame containing study data
#' @param priors Named list of prior specifications
#' @param measure Effect measure type
#' @param section_label Label for this analysis section
#' @param prior_labels Named vector of prior labels
#' @param study_var Study identifier variable (quoted)
#'
#' @return A tibble with posterior draws from mixture models
#'
#' @details
#' Fits a two-component Gaussian mixture model where:
#' \itemize{
#'   \item Component 1: Uses the specified prior
#'   \item Component 2: Uses a wider prior for robustness
#'   \item Weights: Dirichlet(2, 2) prior
#' }
#'
#' @keywords internal
#' @noRd
mixture_model_to_sensitivity_draws <- function(model,
                                               data,
                                               priors,
                                               measure,
                                               section_label,
                                               prior_labels,
                                               study_var) {
  
  study_var_name <- rlang::as_name(study_var)
  
  purrr::imap_dfr(priors, function(prior_obj, prior_name) {
    
    tryCatch({
      
      # -------------------------------------------------
      # 1. Basic checks
      # -------------------------------------------------
      if (!all(c("yi", "sei") %in% names(data))) {
        stop("Mixture model requires columns `yi` and `sei`.")
      }
      
      # -------------------------------------------------
      # 2. Mixture meta-analytic model
      #    (same structure as model_cont)
      # -------------------------------------------------
      study_var_name <- rlang::as_name(study_var)
      
      mix_formula <- brms::bf(
        stats::reformulate(
          termlabels = paste0("(1 | ", study_var_name, ")"),
          response   = "yi | se(sei)"
        ),
        family = brms::mixture(
          stats::gaussian(),
          stats::gaussian()
        )
      )
      
      
      # -------------------------------------------------
      # 3. Extract μ prior from prior list
      # -------------------------------------------------
      mu_prior_vec <- prior_obj |>
        dplyr::filter(class == "Intercept") |>
        dplyr::pull(prior)
      
      mu_prior_str <- mu_prior_vec[1] %||% "normal(0, 1)"
      
      
      # -------------------------------------------------
      # 4. Priors
      #    - μ₁ uses user-supplied prior
      #    - μ₂ is wider (robust component)
      #    - mixture weights use Dirichlet
      #    - τ prior left to brms defaults
      # -------------------------------------------------
      mu_priors <- list(
        brms::set_prior(mu_prior_str, class = "Intercept", dpar = "mu1"),
        brms::set_prior("normal(0, 2)", class = "Intercept", dpar = "mu2")
      )
      
      weight_prior <- brms::set_prior(
        "dirichlet(2, 2)",
        class = "theta"
      )
      
      mix_priors <- do.call(
        c,
        c(mu_priors, list(weight_prior))
      )
      
      # -------------------------------------------------
      # 5. Fit mixture model
      # -------------------------------------------------
      fit <- brms::brm(
        mix_formula,
        data = data,
        prior = mix_priors,
        iter = 4000,
        warmup = 2000,
        chains = 2,
        backend = getOption("brms.backend", "cmdstanr"),
        control = list(adapt_delta = 0.99),
        refresh = 0,
        silent = TRUE
      )
      
      # -------------------------------------------------
      # 6. Mixture-mean posterior (Option A)
      # -------------------------------------------------
      draws <- posterior::as_draws_df(fit)
      
      mu_mix <-
        draws$theta1 * draws$b_mu1_Intercept +
        (1 - draws$theta1) * draws$b_mu2_Intercept
      
      x_draws <-
        if (measure %in% c("OR", "RR", "HR", "IRR")) {
          exp(mu_mix)
        } else {
          mu_mix
        }
      
      tibble::tibble(
        section_label = section_label,
        prior = prior_name,
        prior_label = prior_labels[[prior_name]],
        x = x_draws
      )
      
    }, error = function(e) {
      
      warning(
        sprintf(
          "Failed mixture model for prior %s: %s",
          prior_name, e$message
        ),
        call. = FALSE
      )
      
      tibble::tibble(
        section_label = character(),
        prior = character(),
        prior_label = character(),
        x = numeric()
      )
    })
  })
}

#' PET-PEESE Sensitivity Analysis
#'
#' @description
#' Performs precision-effect test (PET) and precision-effect estimate with
#' standard error (PEESE) for publication bias adjustment. Uses standard errors
#' (SE) as the predictor, not p-values or posterior probabilities.
#'
#' @param model Original brmsfit object
#' @param data Data frame containing study data (must include yi and sei columns)
#' @param priors Named list of prior specifications
#' @param measure Effect measure type (OR, RR, HR, IRR, SMD, MD)
#' @param section_label Label for this analysis section
#' @param prior_labels Named vector of prior labels
#' @param study_var Study identifier variable (quoted)
#' @param direction Direction of expected bias ("negative" or "positive")
#' @param threshold Threshold for choosing between PET and PEESE (default 0.10)
#'
#' @return A tibble with bias-adjusted posterior draws
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Calculates posterior probabilities for each study (to decide PET vs PEESE)
#'   \item Determines whether to use PET (predictor = SE) or PEESE (predictor = SE²)
#'   \item Re-fits model with standard error as bias predictor
#'   \item Returns corrected effect at zero bias (intercept)
#' }
#'
#' PET regression:  effect ~ SE
#' PEESE regression: effect ~ SE²
#' 
#' The bias-corrected estimate is the intercept (effect when SE = 0).
#'
#' @keywords internal
#' @noRd
pet_peese_to_sensitivity_draws <- function(model,
                                           data,
                                           priors,
                                           measure,
                                           section_label,
                                           prior_labels,
                                           study_var,
                                           direction = "negative",
                                           threshold = 0.10) {
  
  # ----------------------------
  # 1. Validation
  # ----------------------------
  if (rlang::quo_is_null(study_var)) {
    stop("`study_var` is required for PET-PEESE analysis", call. = FALSE)
  }
  
  study_var_name <- rlang::as_name(study_var)
  
  if (!study_var_name %in% names(data)) {
    stop(paste("Study variable", study_var_name, "not found in data"), call. = FALSE)
  }
  
  if (!"sei" %in% names(data)) {
    stop("Column 'sei' (standard error) must be present in data for PET-PEESE", call. = FALSE)
  }
  
  if (!"yi" %in% names(data)) {
    stop("Column 'yi' (effect size) must be present in data for PET-PEESE", call. = FALSE)
  }
  
  # ----------------------------
  # 2. Extract study-level posteriors (for PET vs PEESE decision)
  # ----------------------------
  # This is used ONLY to decide whether to use PET or PEESE
  # The actual bias correction uses standard errors, not posteriors
  
  study_posteriors <- extract_study_level_effects(model, data, study_var_name, measure)
  
  # Calculate null value
  null_value <- switch(
    measure,
    OR = 1, RR = 1, HR = 1, IRR = 1,
    MD = 0, SMD = 0
  )
  
  # Calculate posterior probabilities (for decision only)
  data <- data |>
    dplyr::mutate(
      post_prob = purrr::map_dbl(
        seq_len(dplyr::n()),
        ~ {
          if (direction == "negative") {
            mean(study_posteriors[[.x]] < null_value)
          } else {
            mean(study_posteriors[[.x]] > null_value)
          }
        }
      )
    )
  
  # ----------------------------
  # 3. Decide PET or PEESE (using posterior probabilities)
  # ----------------------------
  # If most studies show evidence of an effect, use PEESE
  # Otherwise use PET
  use_peese <- mean(data$post_prob) > threshold
  
  if (!is.na(use_peese)) {
    method_label <- ifelse(use_peese, "PEESE", "PET")
  } else {
    use_peese <- FALSE
    method_label <- "PET"
  }
  
  # ----------------------------
  # 4. Add bias predictor (CORRECTED: using standard errors)
  # ----------------------------
  # KEY CORRECTION: Use standard errors, not posterior probabilities!
  data <- data |>
    dplyr::mutate(
      bias_predictor = if (use_peese) sei^2 else sei
    )
  
  # ----------------------------
  # 5. Fit models with different priors, including bias adjustment
  # ----------------------------
  result <- purrr::imap_dfr(priors, function(prior_obj, prior_name) {
    
    tryCatch({
      
      # ----------------------------
      # 5a. Extract original formula
      # ----------------------------
      original_formula <- formula(model)
      
      # ----------------------------
      # 5b. Add bias predictor to formula
      # ----------------------------
      # The original formula might be:  yi | se(sei) ~ 1 + (1 | study)
      # We need to add: + bias_predictor
      
      formula_bias <- update(
        original_formula,
        . ~ . + bias_predictor
      )
      
      # ----------------------------
      # 5c. Create prior for bias coefficient
      # ----------------------------
      # Prior for the bias coefficient (slope)
      # Weakly informative: allows for positive or negative bias
      prior_bias <- brms::prior(
        "normal(0, 1)",
        class = "b",
        coef = "bias_predictor"
      )
      
      # Combine with user-specified priors
      combined_prior <- c(prior_obj, prior_bias)
      
      # ----------------------------
      # 5d. Fit model with bias correction
      # ----------------------------
      updated_model <- brms::update(
        model,
        formula = formula_bias,
        newdata = data,
        prior = combined_prior,
        refresh = 0,
        silent = 2,
        recompile = FALSE
      )
      
      # ----------------------------
      # 5e. Extract bias-corrected effect
      # ----------------------------
      # The intercept represents the effect when bias_predictor = 0
      # (i.e., when SE = 0, which is the "corrected" estimate)
      
      draws <- tidybayes::spread_draws(updated_model, b_Intercept)
      
      # Transform if necessary (for ratio measures on log scale)
      if (measure %in% c("OR", "RR", "HR", "IRR")) {
        draws$x <- exp(draws$b_Intercept)
      } else {
        draws$x <- draws$b_Intercept
      }
      
      # ----------------------------
      # 5f. Return formatted results
      # ----------------------------
      tibble::tibble(
        section_label = paste0(section_label, " (", method_label, ")"),
        prior = prior_name,
        prior_label = prior_labels[[prior_name]],
        x = draws$x
      )
      
    }, error = function(e) {
      warning(
        sprintf("Failed to fit PET-PEESE for prior %s: %s", prior_name, e$message),
        call. = FALSE
      )
      
      # Return empty tibble with correct structure
      tibble::tibble(
        section_label = character(),
        prior = character(),
        prior_label = character(),
        x = numeric()
      )
    })
  })
  
  # ----------------------------
  # 6. Validation of results
  # ----------------------------
  if (nrow(result) == 0) {
    warning("PET-PEESE produced no results", call. = FALSE)
  }
  
  return(result)
}
