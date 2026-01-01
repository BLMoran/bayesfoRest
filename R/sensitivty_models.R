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
#' standard error (PEESE) for publication bias adjustment.
#'
#' @param model Original brmsfit object
#' @param data Data frame containing study data
#' @param priors Named list of prior specifications
#' @param measure Effect measure type
#' @param section_label Label for this analysis section
#' @param prior_labels Named vector of prior labels
#' @param study_var Study identifier variable (quoted)
#' @param direction Direction of expected bias ("negative" or "positive")
#' @param threshold Threshold for choosing between PET and PEESE
#'
#' @return A tibble with bias-adjusted posterior draws
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Calculates posterior probabilities for each study
#'   \item Determines whether to use PET or PEESE based on threshold
#'   \item Adds bias predictor (SE for PET, SE² for PEESE)
#'   \item Re-fits model with bias adjustment
#'   \item Returns corrected effect at zero bias
#' }
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
  
  # Ensure we have a valid study_var
  if (rlang::quo_is_null(study_var)) {
    stop("study_var is required for PET-PEESE analysis")
  }
  
  # Get the study variable name
  study_var_name <- rlang::as_name(study_var)
  
  # Verify the study variable exists in the data
  if (!study_var_name %in% names(data)) {
    stop(paste("Study variable", study_var_name, "not found in data"))
  }
  
  # Extract study-level posteriors from the original model
  study_posteriors <- extract_study_level_effects(model, data, study_var_name, measure)
  
  # Calculate posterior probabilities for bias indicator
  null_value <- switch(
    measure,
    OR = 1, RR = 1, HR = 1, IRR = 1,
    MD = 0, SMD = 0
  )
  
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
  
  # Determine PET or PEESE
  use_peese <- mean(data$post_prob) > threshold
  
  # Add bias predictor to data
  data <- data |>
    dplyr::mutate(
      bias_predictor = if (use_peese) post_prob^2 else post_prob
    )
  
  # Run sensitivity with different priors, including bias adjustment
  result <- purrr::imap_dfr(priors, function(prior_obj, prior_name) {
    
    tryCatch({
      # Add prior for bias coefficient
      prior_bias <- brms::prior(
        normal(0, 1),
        class = "b",
        coef = "bias_predictor"
      )
      
      combined_prior <- c(prior_obj, prior_bias)
      
      # Update formula to include bias predictor
      formula_bias <- update(
        formula(model),
        . ~ . + bias_predictor
      )
      
      # Fit model with bias correction
      updated_model <- update(
        model,
        formula = formula_bias,
        newdata = data,
        prior = combined_prior,
        refresh = 0,
        recompile = FALSE
      )
      
      # Extract corrected effect (at zero bias)
      draws <- tidybayes::spread_draws(updated_model, b_Intercept)
      
      if (measure %in% c("OR", "RR", "HR", "IRR")) {
        draws$x <- exp(draws$b_Intercept)
      } else {
        draws$x <- draws$b_Intercept
      }
      
      tibble::tibble(
        section_label = paste0(section_label, " (", ifelse(use_peese, "PEESE", "PET"), ")"),
        prior = prior_name,
        prior_label = prior_labels[[prior_name]],
        x = draws$x
      )
    }, error = function(e) {
      warning(paste("Failed to fit PET-PEESE for prior", prior_name, ":", e$message))
      # Return empty tibble with correct structure
      tibble::tibble(
        section_label = character(),
        prior = character(),
        prior_label = character(),
        x = numeric()
      )
    })
  })
  
  # Ensure we always return a tibble
  if (nrow(result) == 0) {
    warning("PET-PEESE produced no results")
  }
  
  return(result)
}
