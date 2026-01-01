#' Internal function for Bayesian Model Averaging for Meta-Analysis using brms
#'
#' @noRd
#' @keywords internal
RoBMA_brms <- function(data, 
                      yi = "yi",
                      sei = "sei",
                      studyvar = NULL,
                      measure = NULL,
                      priors_mu = list(
                        null = "null",
                        normal = brms::prior("normal(0, 1)", class = "Intercept")),
                      priors_tau = list(
                        fixed = "fixed",
                        half_cauchy = brms::prior("cauchy(0, 0.5)", class = "sd")),
                      include_bias = FALSE,
                      use_bridge = TRUE,
                      chains = 4,
                      iter = 4000,
                      warmup = 2000,
                      cores = 4,
                      seed = NULL,
                      silent = TRUE) {
  
  # ----------------------------
  # data prep
  # ----------------------------

  studyvar_quo <- rlang::enquo(studyvar)
  
  if (!rlang::quo_is_null(rlang::enquo(studyvar_quo))) {
    data$study_id <- data[[rlang::as_name(rlang::enquo(studyvar_quo))]]
  }
  
  # ----------------------------
  # model grid
  # ----------------------------
  base_specs <- expand.grid(
    mu   = names(priors_mu),
    tau  = names(priors_tau),
    bias = "none",
    stringsAsFactors = FALSE
  )
  
  if (include_bias) {
    bias_specs <- expand.grid(
      mu   = setdiff(names(priors_mu), "null"),
      tau  = names(priors_tau),
      bias = c("PET", "PEESE"),
      stringsAsFactors = FALSE
    )
    
    model_specs <- rbind(base_specs, bias_specs)
  } else {
    model_specs <- base_specs
  }
  
  # ----------------------------
  # model fitter (MINIMALLY CHANGED)
  # ----------------------------
  fit_one <- function(mu_name, tau_name, bias_name) {
    
    mu_prior  <- priors_mu[[mu_name]]
    tau_prior <- priors_tau[[tau_name]]
    
    is_null <- identical(mu_name, "null")
    
    re_term <- if (!rlang::quo_is_null(studyvar_quo) && tau_name != "fixed") {
      paste0(" + (1 | ", rlang::as_name(studyvar_quo), ")")
    } else ""
    
    bias_term <- switch(
      bias_name,
      PET   = paste0(" + ", sei),
      PEESE = paste0(" + I(", sei, "^2)"),
      ""
    )
    
    formula <- as.formula(
      paste0(
        yi, " | se(", sei, ") ~ ",
        if (is_null) "0 + Intercept" else "1",
        bias_term,
        re_term
      )
    )
    
    prior <- list()
    
    # μ prior
    if (is_null) {
      prior <- append(
        prior,
        list(brms::prior("normal(0, 0.001)", class = "b", coef = "Intercept"))
      )
    } else {
      prior <- append(prior, list(mu_prior))
    }
    
    # τ prior only if random effects
    if (tau_name != "fixed" && !rlang::quo_is_null(studyvar_quo)) {
      prior <- append(prior, list(tau_prior))
    }
    
    # flatten once
    prior <- do.call(c, prior)
    
    
    brms::brm(
      formula = formula,
      data    = data,
      prior   = prior,
      chains  = chains,
      iter    = iter,
      warmup = warmup,
      cores   = cores,
      seed    = seed,
      refresh = 0,
      silent  = silent,
      save_pars = brms::save_pars(all = TRUE),
      control = list(adapt_delta = 0.95)
    )
  }
  
  # ----------------------------
  # fit models
  # ----------------------------
  fits <- vector("list", nrow(model_specs))
  names(fits) <- apply(model_specs, 1, paste, collapse = "_")
  
  for (i in seq_len(nrow(model_specs))) {
    spec <- model_specs[i, ]
    fits[[i]] <- fit_one(spec$mu, spec$tau, spec$bias)
  }
  
  # ----------------------------
  # log marginal likelihoods (THE FIX)
  # ----------------------------
  extract_logml <- function(fit) {
    if (!use_bridge) {
      waic <- tryCatch(brms::waic(fit), error = function(e) NULL)
      if (is.null(waic)) return(NA_real_)
      return(-0.5 * waic$estimates["waic", "Estimate"])
    }
    
    bs <- tryCatch(
      bridgesampling::bridge_sampler(fit, silent = TRUE),
      error = function(e) NULL
    )
    
    if (is.null(bs)) return(NA_real_)
    if (!is.numeric(bs$logml)) return(NA_real_)
    bs$logml
  }
  
  logml <- vapply(fits, extract_logml, numeric(1))
  
  keep <- is.finite(logml)
  
  if (!any(keep)) {
    stop("No models returned valid marginal likelihoods.", call. = FALSE)
  }
  
  fits        <- fits[keep]
  logml       <- logml[keep]
  model_specs <- model_specs[keep, , drop = FALSE]
  
  # ----------------------------
  # weights
  # ----------------------------
  logml <- logml - max(logml)
  weights <- exp(logml) / sum(exp(logml))
  
  # ----------------------------
  # model-averaged posterior
  # ----------------------------
  draws <- purrr::map(fits, ~ as.data.frame(.x)$b_Intercept)
  
  n <- length(draws[[1]])
  ma_draws <- numeric(n)
  
  for (i in seq_len(n)) {
    k <- sample.int(length(draws), 1, prob = weights)
    ma_draws[i] <- draws[[k]][i]
  }
  
  prior_specs_mu <- tibble::tibble(
    mu_prior = names(priors_mu),
    mu_prior_string = purrr::map_chr(
      priors_mu,
      ~ if (is.character(.x) && .x == "null") {
        "normal(0, 0.001)"
      } else {
        as.character(.x$prior)
      }
    )
  ) |>
    dplyr::distinct()
  
  prior_specs_tau <- tibble::tibble(
    tau_prior = names(priors_tau),
    tau_prior_string = purrr::map_chr(
      priors_tau,
      ~ if (is.character(.x) && .x == "fixed") {
        "fixed"
      } else {
        as.character(.x$prior)
      }
    )
  )
  
  
  # ----------------------------
  # return object (RoBMA-compatible)
  # ----------------------------
  out <- list(
    models        = fits,
    model_specs  = model_specs,
    prior_specs = list(
      mu  = prior_specs_mu,
      tau = prior_specs_tau
    ),
    model_weights = weights,
    logml         = logml,
    ma_posterior  = ma_draws,
    measure       = measure
  )
  
  class(out) <- "RoBMA_brms"
  out
}
