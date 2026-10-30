#' Extract Pooled Effect Draws (Backend-Agnostic)
#'
#' @description
#' Extracts posterior draws for the pooled effect (mu) from any supported model
#' backend and transforms them to the outcome scale.
#'
#' @param model A fitted model object (any supported backend).
#' @param data Data frame containing study data.
#' @param measure Effect measure type.
#'
#' @return Numeric vector of posterior draws on the outcome scale.
#'
#' @keywords internal
#' @noRd
extract_pooled_draws <- function(model, data, measure) {

  all_draws <- extract_forest_draws(model, data)

  pooled <- all_draws |>
    dplyr::filter(Author == "Pooled Effect") |>
    dplyr::pull(theta_study)

  # Transform to outcome scale for ratio measures
  if (measure %in% c("OR", "RR", "HR", "IRR")) {
    pooled <- exp(pooled)
  }

  pooled
}

#' Does the fitted model include a random (between-study) effect?
#' @keywords internal
#' @noRd
model_has_re <- function(model) UseMethod("model_has_re")

#' @noRd
model_has_re.default <- function(model) {
  # MetaStan / meta_stan objects carry it in stanDat
  if (!is.null(model$stanDat) && !is.null(model$stanDat$re))
    return(isTRUE(as.logical(model$stanDat$re)))
  NA   # unknown backend
}

#' @noRd
model_has_re.brmsfit <- function(model) {
  any(grepl("^sd_", brms::variables(model)))     # a study-level SD term => random effects
}

#' @noRd
model_has_re.bayesmeta <- function(model) TRUE   # bayesmeta is always random-effects

#' Update Model with New Priors (S3 Generic)
#'
#' @description
#' Refits a model with new prior specifications while keeping the same data
#' and model structure. Dispatches to backend-specific methods. The
#' \code{new_prior} argument is the backend-specific prior object produced by
#' \code{translate_priors_for_backend()}.
#'
#' @param model A fitted model object.
#' @param data Data frame for the model.
#' @param new_prior Backend-specific prior object.
#' @param ... Additional arguments passed to methods.
#'
#' @return Updated model object of the same class.
#'
#' @keywords internal
#' @noRd
update_model_prior <- function(model, data, new_prior, ...) {
  UseMethod("update_model_prior")
}

#' @noRd
update_model_prior.default <- function(model, data, new_prior, ...) {

  rlang::abort(paste0(
    "Prior updating not supported for model class '", class(model)[1], "'.\n",
    "Supported backends: brmsfit, bayesmeta, meta_stan, bayesma, brma, rjags."
  ))
}

#' @noRd
update_model_prior.brmsfit <- function(model, data, new_prior, ...) {
  dots     <- list(...)
  force_re <- dots$force_re

  has_re   <- any(grepl("^sd_", brms::variables(model)))
  want_re  <- if (is.null(force_re)) has_re else isTRUE(force_re)

  form <- model$formula
  if (!want_re && has_re) {
    form <- stats::update(brms::brmsformula(form), . ~ . - (1 | Author))
  }

  prior <- new_prior
  if (!want_re) prior <- prior[prior$class != "sd", , drop = FALSE]

  # Structure change always needs a recompile. cmdstanr fits loaded from disk
  # have no live executable, so never pass recompile = FALSE there -- use NULL
  # (brms/cmdstanr rebuilds only if required).
  structure_changed <- (want_re != has_re)
  backend <- model$backend %||% "rstan"
  recompile <- if (structure_changed) TRUE
  else if (backend == "cmdstanr") NULL     # let it rebuild if needed
  else FALSE                               # rstan can reuse

  stats::update(
    model,
    formula.  = form,
    newdata   = data,
    prior     = prior,
    refresh   = 0,
    silent    = 2,
    recompile = recompile
  )
}

#' @noRd
update_model_prior.bayesmeta <- function(model, data, new_prior, ...) {
  if (!requireNamespace("bayesmeta", quietly = TRUE)) {
    rlang::abort("Package 'bayesmeta' is required.")
  }

  # new_prior is a list with mu_prior (numeric(2)) and tau_prior (function)
  y_col  <- if ("y"  %in% names(data)) "y"  else if ("yi"  %in% names(data)) "yi"  else NULL
  se_col <- if ("se" %in% names(data)) "se" else if ("sei" %in% names(data)) "sei" else NULL

  if (is.null(y_col) || is.null(se_col)) {
    rlang::abort("bayesmeta refit requires 'y'/'yi' and 'se'/'sei' columns.")
  }

  bayesmeta::bayesmeta(
    y         = data[[y_col]],
    sigma     = data[[se_col]],
    labels    = data$Author,
    mu.prior  = new_prior$mu_prior  %||% model$mu.prior,
    tau.prior = new_prior$tau_prior %||% model$tau.prior
  )
}

#' @noRd
update_model_prior.meta_stan <- function(model, data, new_prior, ...) {
  if (!requireNamespace("MetaStan", quietly = TRUE)) {
    rlang::abort("Package 'MetaStan' is required.")
  }

  dots <- list(...)
  force_re <- dots$force_re  # NULL = use original, TRUE/FALSE = override

  # Map link number to likelihood string
  likelihood <- switch(
    as.character(model$stanDat$link),
    "1" = "normal",
    "2" = "binomial",
    "3" = "Poisson",
    "binomial"
  )

  # Determine armVars based on likelihood
  arm_vars <- if (likelihood == "binomial") {
    c(responders = "r", sampleSize = "n")
  } else if (likelihood == "normal") {
    c(mean = "y", sd = "y_se")
  } else {
    c(counts = "count", exposure = "exposure")
  }

  # Filter data_wide to match studies in the (possibly filtered) data
  data_wide <- model$data_wide
  if (!is.null(data) && "Author" %in% names(data) && "Author" %in% names(data_wide)) {
    data_wide <- data_wide[data_wide$Author %in% unique(data$Author), , drop = FALSE]
  }

  stan_dat <- MetaStan::create_MetaStan_dat(
    dat      = data_wide,
    armVars  = arm_vars,
    nArmsVar = "nd"
  )

  # Determine RE flag: force_re overrides, otherwise use original
  re <- if (!is.null(force_re)) force_re else as.logical(model$stanDat$re)

  # Common-effect models (re = FALSE) declare tau_prior as a scalar in Stan;
  # only override tau_prior when the model actually has random effects.
  if (re) {
    tau_prior_val  <- new_prior$tau_prior      %||% model$stanDat$tau_prior
    tau_prior_dist <- new_prior$tau_prior_dist  %||% model$tau_prior_dist
  } else {
    tau_prior_val  <- model$stanDat$tau_prior
    tau_prior_dist <- model$tau_prior_dist
  }

  MetaStan::meta_stan(
    data           = stan_dat,
    likelihood     = likelihood,
    mu_prior       = new_prior$mu_prior  %||% model$stanDat$mu_prior,
    theta_prior    = model$stanDat$theta_prior,
    tau_prior      = tau_prior_val,
    tau_prior_dist = tau_prior_dist,
    beta_prior     = model$stanDat$beta_prior,
    re             = re,
    ncp            = as.logical(model$stanDat$ncp),
    interval.type  = model$interval.type,
    mreg           = as.logical(model$stanDat$mreg)
  )
}


#' Perform Prior Sensitivity Analysis (Backend-Agnostic)
#'
#' @description
#' Re-fits a model with different prior specifications and extracts posterior
#' draws for sensitivity analysis. Priors are specified as \code{bf_prior}
#' objects and translated to the appropriate backend format before refitting.
#'
#' @param model Fitted model object.
#' @param data Data frame containing study data.
#' @param priors Named list of prior specifications, each containing
#'   \code{name}, \code{mu_prior}, and \code{tau_prior}.
#' @param measure Effect measure type.
#' @param section_label Label for this section of the analysis.
#'
#' @return A tibble with columns: section_label, prior, prior_label, x (draws)
#'
#' @keywords internal
#' @noRd
update_prior_sensitivity <- function(model,
                                     data,
                                     priors,
                                     measure,
                                     section_label,
                                     force_re = NULL) {

  re_flag <- if (!is.null(force_re)) isTRUE(force_re) else model_has_re(model)

  purrr::imap_dfr(priors, function(prior_spec, prior_name) {

    prior_label <- prior_spec$name %||% prior_name

    tryCatch({
      backend_prior <- translate_priors_for_backend(prior_spec, model)

      updated_model <- update_model_prior(
        model     = model,
        data      = data,
        new_prior = backend_prior,
        force_re  = force_re
      )

      tibble::tibble(
        section_label = section_label,
        prior         = prior_name,
        prior_label   = prior_label,
        re            = re_flag,
        x             = extract_pooled_draws(updated_model, data, measure)
      )
    }, error = function(e) {
      rlang::warn(paste0(
        "Failed prior sensitivity for '", prior_name, "' in '",
        section_label, "': ", e$message
      ))
      tibble::tibble(
        section_label = character(),
        prior         = character(),
        prior_label   = character(),
        re            = logical(),
        x             = numeric()
      )
    })
  })
}


#' Summarise Sensitivity Posteriors
#'
#' @description
#' Computes summary statistics and probabilities from posterior draws across
#' all sensitivity analyses.
#'
#' @param draws Tibble containing posterior draws.
#' @param null_value Null value for hypothesis testing.
#' @param null_range Optional range for region of practical equivalence.
#'
#' @return Tibble with summary statistics.
#'
#' @keywords internal
#' @noRd
summarise_sensitivity_posteriors <- function(draws, null_value, null_range) {

  section_levels <- unique(draws$section_label)
  prior_levels   <- unique(draws$prior)

  draws |>
    dplyr::mutate(
      section_label = factor(section_label, levels = section_levels),
      prior         = factor(prior, levels = prior_levels)
    ) |>
    dplyr::group_by(section_label, prior, prior_label, re) |>
    dplyr::summarise(
      median = stats::median(x),
      l95    = stats::quantile(x, 0.025),
      u95    = stats::quantile(x, 0.975),

      pr_benefit    = round(mean(x < null_value) * 100, 1),
      pr_no_benefit = round(mean(x > null_value) * 100, 1),

      pr_benefit_null_range =
        if (!is.null(null_range)) round(mean(x < null_range[1]) * 100, 1) else NA_real_,
      pr_no_benefit_null_range =
        if (!is.null(null_range)) round(mean(x > null_range[2]) * 100, 1) else NA_real_,

      .groups = "drop"
    ) |>
    dplyr::arrange(section_label, prior)
}


#' Run MA Sensitivity Dispatcher
#'
#' @description
#' Dispatcher function that runs the appropriate sensitivity analysis based on
#' the specification provided.
#'
#' @param spec List containing analysis specifications.
#' @param priors Named list of prior specifications.
#' @param measure Effect measure type.
#'
#' @return Tibble with posterior draws and labels.
#'
#' @keywords internal
#' @noRd
run_ma_sensitivity <- function(spec, priors, measure) {

  update_prior_sensitivity(
    model         = spec$model,
    data          = spec$data,
    priors        = priors,
    measure       = measure,
    section_label = spec$section_label,
    force_re      = spec$force_re  # NULL unless explicitly set
  )
}
