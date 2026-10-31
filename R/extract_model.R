#' Model Extraction Functions for bayesfoRest
#'
#' @description
#' S3 generics and methods for extracting posterior draws, pooled effect
#' summaries, and model type information from various Bayesian meta-analysis
#' backends. Supported backends: brms, MetaStan, bayesmeta, rstanarm, rstan,
#' cmdstanr, bayesma, RoBMA (brma), R2jags (JAGS).
#'
#' All methods produce output in a canonical internal format so that
#' downstream plotting and table functions are backend-agnostic.
#'
#' @name extract-model
#' @keywords internal
NULL


# =============================================================================
# CANONICAL DRAWS FORMAT
# =============================================================================
#
# Every extract_forest_draws() method must return a tibble with columns:
#
#   .draw        : integer draw index
#   Author       : character study name (includes "Pooled Effect" rows)
#   theta_study  : numeric study-level effect on model scale
#   eta          : numeric random deviation from pooled (NA for pooled / CE)
#   tau          : numeric between-study heterogeneity SD (NA for CE models)
#


# =============================================================================
# detect_model_info
# =============================================================================

#' Detect Model Backend and Pooling Type
#'
#' @description
#' Inspects a fitted model object and returns its backend type and whether it
#' uses common-effects or random-effects pooling.
#'
#' @param model A fitted model object.
#' @param ... Additional arguments passed to methods.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{backend}{Character string identifying the backend.}
#'     \item{pooling}{Character string: "random_effect" or "common_effect".}
#'   }
#'
#' @keywords internal
#' @export
detect_model_info <- function(model, ...) {
  UseMethod("detect_model_info")
}

#' @export
detect_model_info.default <- function(model, ...) {
  rlang::abort(
    c(
      paste0("Unsupported model class: '", class(model)[1], "'"),
      "i" = "Supported types: brmsfit, meta_stan, bayesmeta, stanreg, stanfit, CmdStanMCMC, bayesma, brma, rjags"
    )
  )
}

#' @export
detect_model_info.brmsfit <- function(model, ...) {
  posterior_df <- as.data.frame(model)
  has_re <- any(grepl("^sd_", names(posterior_df)))
  list(backend = "brms",
       pooling = if (has_re) "random_effect" else "common_effect")
}

#' @export
detect_model_info.meta_stan <- function(model, ...) {
  posterior_df <- as.data.frame(model$fit)
  has_tau <- any(grepl("^tau", names(posterior_df)))
  list(backend = "metastan",
       pooling = if (has_tau) "random_effect" else "common_effect")
}

#' @export
detect_model_info.bayesmeta <- function(model, ...) {
  list(backend = "bayesmeta", pooling = "random_effect")
}

#' @export
detect_model_info.stanreg <- function(model, ...) {
  posterior_df <- as.data.frame(model)
  has_re <- any(grepl("^b\\[|^Sigma\\[|^sd_", names(posterior_df)))
  list(backend = "rstanarm",
       pooling = if (has_re) "random_effect" else "common_effect")
}

#' @export
detect_model_info.stanfit <- function(model, tau_par = "tau", ...) {
  posterior_df <- as.data.frame(model)
  has_tau <- any(grepl(paste0("^", tau_par, "$|^", tau_par, "\\["),
                       names(posterior_df)))
  list(backend = "rstan",
       pooling = if (has_tau) "random_effect" else "common_effect")
}

#' @export
detect_model_info.CmdStanMCMC <- function(model, tau_par = "tau", ...) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    rlang::abort("Package 'posterior' is required for CmdStanMCMC models.")
  }
  posterior_df <- as.data.frame(posterior::as_draws_df(model$draws()))
  has_tau <- any(grepl(paste0("^", tau_par, "$|^", tau_par, "\\["),
                       names(posterior_df)))
  list(backend = "cmdstanr",
       pooling = if (has_tau) "random_effect" else "common_effect")
}


# =============================================================================
# extract_pooled_summary
# =============================================================================

#' Extract Pooled Effect Summary
#'
#' @description
#' Extracts the pooled (overall) effect estimate and credible interval from a
#' fitted model. Values are on the model scale (log scale for ratio measures).
#'
#' @param model A fitted model object.
#' @param prob Probability mass for the credible interval (default 0.95).
#' @param ... Additional arguments passed to methods.
#'
#' @return A named numeric vector: \code{estimate}, \code{lower}, \code{upper}.
#'
#' @keywords internal
#' @export
extract_pooled_summary <- function(model, prob = 0.95, ...) {
  UseMethod("extract_pooled_summary")
}

#' @export
extract_pooled_summary.default <- function(model, prob = 0.95, ...) {
  rlang::abort(
    c(
      paste0("No extract_pooled_summary method for class '", class(model)[1], "'"),
      "i" = "Supported types: brmsfit, meta_stan, bayesmeta, stanreg, stanfit, CmdStanMCMC, bayesma, brma, rjags"
    )
  )
}

#' @export
extract_pooled_summary.brmsfit <- function(model, prob = 0.95, ...) {
  posterior_df <- as.data.frame(model)
  theta <- posterior_df[["b_Intercept"]]
  if (is.null(theta)) {
    rlang::abort("Could not find b_Intercept parameter in brms model.")
  }
  alpha <- (1 - prob) / 2
  c(estimate = stats::median(theta),
    lower    = unname(stats::quantile(theta, alpha)),
    upper    = unname(stats::quantile(theta, 1 - alpha)))
}

#' @export
extract_pooled_summary.meta_stan <- function(model, prob = 0.95, ...) {
  theta <- as.data.frame(model$fit)[["theta"]]
  alpha <- (1 - prob) / 2
  c(estimate = stats::median(theta),
    lower    = unname(stats::quantile(theta, alpha)),
    upper    = unname(stats::quantile(theta, 1 - alpha)))
}

#' @export
extract_pooled_summary.bayesmeta <- function(model, prob = 0.95, ...) {
  # Use qposterior to get quantiles of mu directly (no MCMC draws)
  alpha <- (1 - prob) / 2
  c(estimate = model$qposterior(mu.p = 0.5),
    lower    = model$qposterior(mu.p = alpha),
    upper    = model$qposterior(mu.p = 1 - alpha))
}

#' @export
extract_pooled_summary.stanreg <- function(model, prob = 0.95, ...) {
  posterior_df <- as.data.frame(model)
  int_col <- grep("^\\(Intercept\\)$|^Intercept$", names(posterior_df), value = TRUE)
  if (length(int_col) == 0) {
    rlang::abort("Could not find intercept parameter in rstanarm model.")
  }
  theta <- posterior_df[[int_col[1]]]
  alpha <- (1 - prob) / 2
  c(estimate = stats::median(theta),
    lower    = unname(stats::quantile(theta, alpha)),
    upper    = unname(stats::quantile(theta, 1 - alpha)))
}

#' @export
extract_pooled_summary.stanfit <- function(model, prob = 0.95,
                                           theta_par = "theta", ...) {
  posterior_df <- as.data.frame(model)
  theta_col <- .find_parameter(posterior_df, theta_par, c("mu", "beta", "d"))
  theta <- posterior_df[[theta_col]]
  alpha <- (1 - prob) / 2
  c(estimate = stats::median(theta),
    lower    = unname(stats::quantile(theta, alpha)),
    upper    = unname(stats::quantile(theta, 1 - alpha)))
}

#' @export
extract_pooled_summary.CmdStanMCMC <- function(model, prob = 0.95,
                                                theta_par = "theta", ...) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    rlang::abort("Package 'posterior' is required for CmdStanMCMC models.")
  }
  posterior_df <- as.data.frame(posterior::as_draws_df(model$draws()))
  theta_col <- .find_parameter(posterior_df, theta_par, c("mu", "beta", "d"))
  theta <- posterior_df[[theta_col]]
  alpha <- (1 - prob) / 2
  c(estimate = stats::median(theta),
    lower    = unname(stats::quantile(theta, alpha)),
    upper    = unname(stats::quantile(theta, 1 - alpha)))
}


# =============================================================================
# extract_forest_draws
# =============================================================================

#' Extract Posterior Draws in Canonical Format
#'
#' @description
#' Extracts study-level and pooled posterior draws from a fitted model,
#' returning them in the canonical format expected by bayesfoRest plotting
#' and table functions.
#'
#' @param model A fitted model object.
#' @param data Data frame with study information. Must include an
#'   \code{Author} column.
#' @param ... Additional arguments passed to methods (e.g., \code{theta_par},
#'   \code{tau_par}, \code{study_par} for stanfit/CmdStanMCMC models).
#'
#' @return A tibble with columns: \code{.draw}, \code{Author},
#'   \code{theta_study}, \code{eta}, \code{tau}.
#'
#' @keywords internal
#' @export
extract_forest_draws <- function(model, data, ...) {
  UseMethod("extract_forest_draws")
}

#' @export
extract_forest_draws.default <- function(model, data, ...) {
  rlang::abort(
    c(
      paste0("No extract_forest_draws method for class '", class(model)[1], "'"),
      "i" = "Supported types: brmsfit, meta_stan, bayesmeta, stanreg, stanfit, CmdStanMCMC, bayesma, brma, rjags"
    )
  )
}

#' @export
extract_forest_draws.brmsfit <- function(model, data, ...) {
  if (!"Author" %in% names(data)) {
    rlang::abort("Data must contain an 'Author' column.")
  }

  n_studies    <- nrow(data)
  study_names  <- as.character(data$Author)
  posterior_df <- as.data.frame(model)
  n_draws      <- nrow(posterior_df)
  draw_ids     <- seq_len(n_draws)

  intercept <- posterior_df[["b_Intercept"]]
  if (is.null(intercept)) {
    rlang::abort("Could not find b_Intercept parameter in brms model.")
  }

  re_cols <- grep("^r_[^\\[]+\\[", names(posterior_df), value = TRUE)

  tau_col <- grep("^sd_[^_]+__Intercept$", names(posterior_df), value = TRUE)
  tau_vec <- if (length(tau_col) > 0) posterior_df[[tau_col[1]]] else rep(NA_real_, n_draws)

  study_draws <- purrr::map_dfr(seq_len(n_studies), function(i) {
    name_escaped <- gsub("([.()\\[\\]])", "\\\\\\1", study_names[i])
    name_brms    <- gsub(" ", ".", name_escaped)
    re_col <- grep(paste0("\\[", name_brms, ","), re_cols, value = TRUE)
    deviation <- if (length(re_col) > 0) posterior_df[[re_col[1]]] else rep(0, n_draws)

    tibble::tibble(
      .draw       = draw_ids,
      Author      = study_names[i],
      theta_study = intercept + deviation,
      eta         = deviation,
      tau         = tau_vec
    )
  })

  pooled_draws <- tibble::tibble(
    .draw       = draw_ids,
    Author      = "Pooled Effect",
    theta_study = intercept,
    eta         = NA_real_,
    tau         = tau_vec
  )

  dplyr::bind_rows(study_draws, pooled_draws)
}

#' @export
extract_forest_draws.meta_stan <- function(model, data, ...) {
  if (!"Author" %in% names(data)) {
    rlang::abort("Data must contain an 'Author' column.")
  }

  posterior_df <- as.data.frame(model$fit)
  n_draws     <- nrow(posterior_df)
  n_studies   <- nrow(data)
  study_names <- as.character(data$Author)
  theta       <- posterior_df[["theta"]]
  draw_ids    <- seq_len(n_draws)

  tau_col <- grep("^tau\\[1\\]$|^tau$", names(posterior_df), value = TRUE)
  tau_vec <- if (length(tau_col) > 0) posterior_df[[tau_col[1]]] else rep(NA_real_, n_draws)

  has_u <- any(grepl("^u\\[1,\\d+\\]$", names(posterior_df)))
  has_d <- any(grepl("^d\\[\\d+\\]$", names(posterior_df)))

  # Detect arm-level vs study-level d indexing
  # MetaStan uses arm-level: 2 arms per study → n_d = 2 * n_studies
  # Control arms at odd indices d[1], d[3], ...; treatment at even d[2], d[4], ...
  if (has_d) {
    d_cols <- grep("^d\\[\\d+\\]$", names(posterior_df), value = TRUE)
    n_d    <- length(d_cols)
    arm_level_d <- (n_d >= 2 * n_studies)
  }

  study_draws <- purrr::map_dfr(seq_len(n_studies), function(i) {
    if (has_u) {
      # Prefer u parameters: study-indexed, non-centered (u ~ N(0,1))
      # Study-specific treatment effect = theta + u[1,i] * tau
      u_col <- paste0("u[1,", i, "]")
      if (u_col %in% names(posterior_df)) {
        deviation    <- posterior_df[[u_col]] * tau_vec
        study_effect <- theta + deviation
      } else {
        deviation    <- rep(0, n_draws)
        study_effect <- theta
      }
    } else if (has_d) {
      if (arm_level_d) {
        # Arm-level indexing: treatment effect = d[treatment] - d[control]
        # For study i: control arm = d[2*i - 1], treatment arm = d[2*i]
        d_ctrl <- paste0("d[", 2 * i - 1, "]")
        d_trt  <- paste0("d[", 2 * i, "]")
        if (d_trt %in% names(posterior_df) && d_ctrl %in% names(posterior_df)) {
          study_effect <- posterior_df[[d_trt]] - posterior_df[[d_ctrl]]
          deviation    <- study_effect - theta
        } else {
          study_effect <- theta
          deviation    <- rep(0, n_draws)
        }
      } else {
        # Study-level indexing (rare, but handle it)
        d_col <- paste0("d[", i, "]")
        if (d_col %in% names(posterior_df)) {
          study_effect <- posterior_df[[d_col]]
          deviation    <- study_effect - theta
        } else {
          study_effect <- theta
          deviation    <- rep(0, n_draws)
        }
      }
    } else {
      # No study-level parameters: common-effects, use theta for all
      deviation    <- rep(0, n_draws)
      study_effect <- theta
    }

    tibble::tibble(
      .draw       = draw_ids,
      Author      = study_names[i],
      theta_study = study_effect,
      eta         = deviation,
      tau         = tau_vec
    )
  })

  pooled_draws <- tibble::tibble(
    .draw       = draw_ids,
    Author      = "Pooled Effect",
    theta_study = theta,
    eta         = NA_real_,
    tau         = tau_vec
  )

  dplyr::bind_rows(study_draws, pooled_draws)
}

#' @export
extract_forest_draws.bayesmeta <- function(model, data, ...) {
  if (!"Author" %in% names(data)) {
    rlang::abort("Data must contain an 'Author' column.")
  }

  n_draws     <- 4000L
  n_studies   <- model$k
  study_names <- as.character(data$Author)
  draw_ids    <- seq_len(n_draws)

  # Fast grid-based sampling of (mu, tau) from the joint posterior.
  # bayesmeta's rposterior() uses inverse-CDF sampling (uniroot + integrate

  # per draw) which is extremely slow. Instead:
  #   1. Sample tau from its marginal posterior via a fine grid
  #   2. Sample mu | tau analytically (conjugate Normal)
  y_obs     <- model$y
  sigma_obs <- model$sigma

  # Build grid for marginal posterior of tau
  tau_upper <- model$qposterior(tau.p = 0.9999)
  tau_grid  <- seq(0, tau_upper, length.out = 2000L)
  tau_dens  <- model$dposterior(tau = tau_grid)
  tau_dens[is.na(tau_dens) | tau_dens < 0] <- 0

  # Sample tau from grid via inverse-CDF (vectorised, instant)
  tau_cdf   <- cumsum(tau_dens)
  tau_cdf   <- tau_cdf / tau_cdf[length(tau_cdf)]
  u         <- stats::runif(n_draws)
  tau_idx   <- findInterval(u, tau_cdf) + 1L
  tau_idx   <- pmin(tau_idx, length(tau_grid))
  tau_draws <- tau_grid[tau_idx]

  # Sample mu | tau analytically: mu | tau, y ~ N(mu_hat, 1/prec_mu)
  prec_mu <- vapply(tau_draws, function(t) {
    sum(1 / (sigma_obs^2 + t^2))
  }, numeric(1))
  mu_hat <- vapply(seq_along(tau_draws), function(j) {
    sum(y_obs / (sigma_obs^2 + tau_draws[j]^2)) / prec_mu[j]
  }, numeric(1))
  mu_draws <- stats::rnorm(n_draws, mean = mu_hat, sd = sqrt(1 / prec_mu))

  # Compute study-specific theta_i analytically from posterior conditional:
  #   theta_i | y_i, sigma_i, mu, tau ~ N(w_i, v_i)
  #   v_i = 1 / (1/sigma_i^2 + 1/tau^2)
  #   w_i = v_i * (y_i/sigma_i^2 + mu/tau^2)
  study_draws <- purrr::map_dfr(seq_len(n_studies), function(i) {
    prec_data  <- 1 / sigma_obs[i]^2
    # When tau = 0, the study effect equals the pooled effect (no heterogeneity)
    tau_safe   <- pmax(tau_draws, .Machine$double.eps)
    prec_prior <- 1 / tau_safe^2
    v_i        <- 1 / (prec_data + prec_prior)
    w_i        <- v_i * (y_obs[i] * prec_data + mu_draws * prec_prior)
    theta_i    <- stats::rnorm(n_draws, mean = w_i, sd = sqrt(v_i))
    # For tau = 0 draws, theta_i should equal mu (complete shrinkage)
    theta_i[tau_draws == 0] <- mu_draws[tau_draws == 0]
    tibble::tibble(
      .draw       = draw_ids,
      Author      = study_names[i],
      theta_study = theta_i,
      eta         = theta_i - mu_draws,
      tau         = tau_draws
    )
  })

  pooled_draws <- tibble::tibble(
    .draw       = draw_ids,
    Author      = "Pooled Effect",
    theta_study = mu_draws,
    eta         = NA_real_,
    tau         = tau_draws
  )

  dplyr::bind_rows(study_draws, pooled_draws)
}

#' @export
extract_forest_draws.stanreg <- function(model, data, ...) {
  if (!"Author" %in% names(data)) {
    rlang::abort("Data must contain an 'Author' column.")
  }

  n_studies    <- nrow(data)
  study_names  <- as.character(data$Author)
  posterior_df <- as.data.frame(model)
  n_draws      <- nrow(posterior_df)
  draw_ids     <- seq_len(n_draws)

  int_col <- grep("^\\(Intercept\\)$|^Intercept$", names(posterior_df), value = TRUE)
  if (length(int_col) == 0) {
    rlang::abort("Could not find intercept parameter in rstanarm model.")
  }
  intercept <- posterior_df[[int_col[1]]]

  re_cols <- grep("^b\\[\\(Intercept\\) Author:", names(posterior_df), value = TRUE)

  tau_col <- grep("^Sigma\\[Author:\\(Intercept\\),\\(Intercept\\)\\]",
                  names(posterior_df), value = TRUE)
  if (length(tau_col) > 0) {
    tau_vec <- sqrt(posterior_df[[tau_col[1]]])
  } else {
    tau_col <- grep("sd_Author", names(posterior_df), value = TRUE)
    tau_vec <- if (length(tau_col) > 0) posterior_df[[tau_col[1]]] else rep(NA_real_, n_draws)
  }

  study_draws <- purrr::map_dfr(seq_len(n_studies), function(i) {
    name_pattern <- gsub(" ", "\\\\.", gsub("\\.", "\\\\.", study_names[i]))
    re_col <- grep(paste0("Author:", name_pattern, "\\]"), re_cols, value = TRUE)
    deviation <- if (length(re_col) > 0) posterior_df[[re_col[1]]] else rep(0, n_draws)

    tibble::tibble(
      .draw       = draw_ids,
      Author      = study_names[i],
      theta_study = intercept + deviation,
      eta         = deviation,
      tau         = tau_vec
    )
  })

  pooled_draws <- tibble::tibble(
    .draw       = draw_ids,
    Author      = "Pooled Effect",
    theta_study = intercept,
    eta         = NA_real_,
    tau         = tau_vec
  )

  dplyr::bind_rows(study_draws, pooled_draws)
}

#' @export
extract_forest_draws.stanfit <- function(model, data,
                                         theta_par = "theta",
                                         tau_par = "tau",
                                         study_par = NULL, ...) {
  if (!"Author" %in% names(data)) {
    rlang::abort("Data must contain an 'Author' column.")
  }
  posterior_df <- as.data.frame(model)
  .extract_stan_draws(posterior_df, data, theta_par, tau_par, study_par)
}

#' @export
extract_forest_draws.CmdStanMCMC <- function(model, data,
                                              theta_par = "theta",
                                              tau_par = "tau",
                                              study_par = NULL, ...) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    rlang::abort("Package 'posterior' is required for CmdStanMCMC models.")
  }
  if (!"Author" %in% names(data)) {
    rlang::abort("Data must contain an 'Author' column.")
  }
  posterior_df <- as.data.frame(posterior::as_draws_df(model$draws()))
  .extract_stan_draws(posterior_df, data, theta_par, tau_par, study_par)
}


# =============================================================================
# REFIT FOR SUBGROUPS
# =============================================================================

#' Refit Model for Subgroup Data
#'
#' Attempts to refit a model on a subset of data for subgroup analysis.
#' Returns the refitted model on success, or NULL if the backend does not
#' support automatic refitting.
#'
#' @param model A fitted model object.
#' @param data Subset of data for this subgroup.
#' @param ... Additional arguments.
#'
#' @return A refitted model, or NULL.
#'
#' @keywords internal
#' @noRd
refit_for_subgroup <- function(model, data, ...) {
  UseMethod("refit_for_subgroup")
}

#' @export
refit_for_subgroup.default <- function(model, data, ...) {
  rlang::warn(
    c(
      paste0("Cannot automatically refit '", class(model)[1],
             "' models for subgroup analysis."),
      "i" = "Subgroup pooled effects will be approximated from the overall model."
    ),
    .frequency = "once",
    .frequency_id = paste0("refit_subgroup_", class(model)[1])
  )
  NULL
}

#' @export
refit_for_subgroup.brmsfit <- function(model, data, ...) {
  md   <- model$data                              # yi, sei, Author
  keep <- as.character(data$Author)               # the studies IN THIS SUBGROUP (already split by age_group upstream)
  sub  <- md[as.character(md$Author) %in% keep, , drop = FALSE]
  if (nrow(sub) < 1) return(NULL)
  tryCatch(
    stats::update(model, newdata = sub, recompile = FALSE, refresh = 0),
    error = function(e) { rlang::warn(conditionMessage(e)); NULL }
  )
}

#' @export
refit_for_subgroup.bayesmeta <- function(model, data, ...) {
  if (!requireNamespace("bayesmeta", quietly = TRUE)) {
    rlang::abort("Package 'bayesmeta' is required.")
  }
  y_col  <- if ("y" %in% names(data)) "y" else if ("yi" %in% names(data)) "yi" else NULL
  se_col <- if ("se" %in% names(data)) "se" else if ("sei" %in% names(data)) "sei" else NULL
  if (is.null(y_col) || is.null(se_col)) {
    rlang::warn("Cannot refit bayesmeta for subgroup: missing y/se columns.")
    return(NULL)
  }
  # Use the original model's tau prior for subgroup refits. When the
  # original prior is proper this preserves the user's specification.
  # When it is improper (the bayesmeta default), few studies per
  # subgroup leave tau unidentified, producing extreme posterior draws
  # (tau > 1000) that make density slabs invisible. In that case we
  # fall back to a vague but proper half-normal(0, 1) prior.
  if (isTRUE(model$tau.prior.proper)) {
    tau_prior <- model$tau.prior
  } else {
    tau_prior <- function(t) bayesmeta::dhalfnormal(t, scale = 1)
  }
  tryCatch(
    bayesmeta::bayesmeta(y = data[[y_col]], sigma = data[[se_col]],
                         labels = data$Author,
                         tau.prior = tau_prior,
                         mu.prior  = model$mu.prior),
    error = function(e) {
      rlang::warn(paste0("bayesmeta subgroup refit failed: ", e$message))
      NULL
    }
  )
}

#' @export
refit_for_subgroup.meta_stan <- function(model, data, ...) {
  if (!requireNamespace("MetaStan", quietly = TRUE)) {
    rlang::abort("Package 'MetaStan' is required.")
  }

  # Map study names to row indices in the original data_wide
  # data_wide$study is 1:n_studies, matching the row order of the user's data
  # The 'data' argument here is the subgroup subset of the user's data frame
  sg_authors <- as.character(data$Author)

  # We need the full original data to find which row indices the subgroup studies are
  # The model stores data_wide with numeric study IDs matching the original data order
  # We need to identify which rows of data_wide correspond to the subgroup
  # The caller passes the subgroup subset, so we match by the row indices
  # stored in the ... or by matching Author names to the original data order

  # Get original study count and subgroup indices
  n_orig <- nrow(model$data_wide)
  n_sg   <- nrow(data)

  # The subgroup data rows must map to data_wide rows

  # We rely on the original_data passed via ... to find indices
  original_data <- list(...)$original_data
  if (!is.null(original_data)) {
    all_authors <- as.character(original_data$Author)
    sg_indices  <- which(all_authors %in% sg_authors)
  } else {
    # Fallback: assume data rows are in same order as data_wide
    # This works when forest_data.R passes the nested subgroup data
    # which preserves row ordering from the original data
    sg_indices <- seq_len(n_sg)
    rlang::warn(
      "Subgroup study indices inferred by position. Pass original_data for exact matching.",
      .frequency = "once",
      .frequency_id = "metastan_subgroup_index"
    )
  }

  # Filter data_wide to subgroup rows and renumber
  sg_data_wide <- model$data_wide[sg_indices, , drop = FALSE]
  sg_data_wide$study <- seq_len(n_sg)

  # Map link number to likelihood string
  likelihood <- switch(
    as.character(model$stanDat$link),
    "1" = "normal",
    "2" = "binomial",
    "3" = "Poisson",
    "binomial"  # fallback
  )

  # Determine armVars based on likelihood
  arm_vars <- if (likelihood == "binomial") {
    c(responders = "r", sampleSize = "n")
  } else if (likelihood == "normal") {
    c(mean = "y", sd = "y_se")
  } else {
    c(counts = "count", exposure = "exposure")
  }

  tryCatch({
    sg_stan_dat <- MetaStan::create_MetaStan_dat(
      dat = sg_data_wide,
      armVars = arm_vars,
      nArmsVar = "nd"
    )

    MetaStan::meta_stan(
      data         = sg_stan_dat,
      likelihood   = likelihood,
      mu_prior     = model$stanDat$mu_prior,
      theta_prior  = model$stanDat$theta_prior,
      tau_prior    = model$stanDat$tau_prior,
      tau_prior_dist = model$tau_prior_dist,
      beta_prior   = model$stanDat$beta_prior,
      re           = as.logical(model$stanDat$re),
      ncp          = as.logical(model$stanDat$ncp),
      interval.type = model$interval.type,
      mreg         = as.logical(model$stanDat$mreg)
    )
  }, error = function(e) {
    rlang::warn(paste0("MetaStan subgroup refit failed: ", e$message))
    NULL
  })
}

#' @export
refit_for_subgroup.stanreg <- function(model, data, ...) {
  tryCatch(
    stats::update(model, data = data, refresh = 0),
    error = function(e) {
      rlang::warn(paste0("rstanarm subgroup refit failed: ", e$message))
      NULL
    }
  )
}


#' Add prediction draws to a canonical draws tibble
#'
#' Generates a predictive distribution for a hypothetical new study by
#' sampling theta_pred = theta_pooled + N(0, tau) for each posterior draw.
#' Works with any backend since it operates on the canonical format.
#'
#' @param draws Canonical draws tibble (output of extract_forest_draws).
#' @return The draws tibble with an added "Prediction" Author.
#' @noRd
add_prediction_draws <- function(draws) {
  pooled <- draws |> dplyr::filter(Author == "Pooled Effect")

  # Need tau for prediction; skip if CE model (all tau = NA)
  if (all(is.na(pooled$tau))) {
    rlang::warn("Cannot add prediction interval for common-effects models (no tau).")
    return(draws)
  }

  set.seed(12345)
  pred_draws <- pooled |>
    dplyr::transmute(
      .draw       = .draw,
      Author      = "Prediction",
      theta_study = theta_study + stats::rnorm(dplyr::n(), mean = 0, sd = tau),
      eta         = NA_real_,
      tau         = tau
    )

  dplyr::bind_rows(draws, pred_draws)
}


# =============================================================================
# bayesma METHODS
# =============================================================================

#' @export
detect_model_info.bayesma <- function(model, ...) {
  has_re <- model$meta$model_type != "common_effect"
  list(backend = "bayesma",
       pooling = if (has_re) "random_effect" else "common_effect")
}

#' @export
extract_pooled_summary.bayesma <- function(model, prob = 0.95, ...) {
  mu_draws <- as.numeric(model$draws$mu$mu)
  alpha <- (1 - prob) / 2
  c(estimate = stats::median(mu_draws),
    lower    = unname(stats::quantile(mu_draws, alpha)),
    upper    = unname(stats::quantile(mu_draws, 1 - alpha)))
}

#' @export
extract_forest_draws.bayesma <- function(model, data, ...) {
  if (!"Author" %in% names(data)) {
    rlang::abort("Data must contain an 'Author' column.")
  }

  is_re    <- model$meta$model_type != "common_effect"
  mu_draws <- as.numeric(model$draws$mu$mu)
  n_draws  <- length(mu_draws)
  draw_ids <- seq_len(n_draws)

  if (is_re && !is.null(model$draws$study)) {
    tau_draws <- as.numeric(model$draws$tau$tau)
    study_col <- model$meta$studyvar %||% "Author"

    study_draws <- model$draws$study |>
      dplyr::transmute(
        .draw       = .draw,
        Author      = as.character(.data[[study_col]]),
        theta_study = theta,
        eta         = theta - mu,
        tau         = tau
      )
  } else {
    tau_draws   <- rep(NA_real_, n_draws)
    study_names <- as.character(data$Author)

    study_draws <- purrr::map_dfr(study_names, function(nm) {
      tibble::tibble(
        .draw       = draw_ids,
        Author      = nm,
        theta_study = mu_draws,
        eta         = 0,
        tau         = tau_draws
      )
    })
  }

  pooled_draws <- tibble::tibble(
    .draw       = draw_ids,
    Author      = "Pooled Effect",
    theta_study = mu_draws,
    eta         = NA_real_,
    tau         = tau_draws
  )

  dplyr::bind_rows(study_draws, pooled_draws)
}

#' @export
refit_for_subgroup.bayesma <- function(model, data, ..., re_min_k = NULL) {
  # bayesma isn't a CRAN dependency, so the package name is kept out of a
  # literal requireNamespace()/"::" call -- both are flagged by R CMD check
  # as undeclared dependencies even when guarded like this.
  pkg <- "bayesma"
  if (!requireNamespace(pkg, quietly = TRUE)) {
    rlang::abort(paste0("Package '", pkg, "' is required."))
  }
  if (is.null(model$call)) {
    rlang::warn("Cannot refit bayesma: no stored call.")
    return(NULL)
  }

  tryCatch({
    # bayes_forest() renames data columns to its own canonical names; point
    # the stored call's column arguments (tidy-eval symbols) at those names
    forest_col_map <- list(
      studyvar   = quote(Author),
      yearvar    = quote(Year),
      n_ctrl     = quote(N_Control),
      n_int      = quote(N_Intervention),
      mean_ctrl  = quote(Mean_Control),
      mean_int   = quote(Mean_Intervention),
      sd_ctrl    = quote(SD_Control),
      sd_int     = quote(SD_Intervention),
      event_ctrl = quote(Event_Control),
      event_int  = quote(Event_Intervention)
    )
    call_args <- as.list(model$call)[-1]
    present <- names(forest_col_map)[
      names(forest_col_map) %in% names(call_args) &
        purrr::map_lgl(forest_col_map, ~ as.character(.x) %in% names(data))
    ]
    new_call <- rlang::call_modify(model$call, !!!forest_col_map[present], data = data)
    eval(new_call, asNamespace(pkg))
  }, error = function(e) {
    rlang::warn(paste0("bayesma subgroup refit failed: ", e$message))
    NULL
  })
}


# =============================================================================
# brma (RoBMA) METHODS
# =============================================================================

#' @export
detect_model_info.brma <- function(model, ...) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    rlang::abort("Package 'posterior' is required for brma models.")
  }
  draws_df <- posterior::as_draws_df(model$fit$mcmc)
  has_tau <- "tau" %in% names(draws_df)
  list(backend = "brma",
       pooling = if (has_tau) "random_effect" else "common_effect")
}

#' @export
extract_pooled_summary.brma <- function(model, prob = 0.95, ...) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    rlang::abort("Package 'posterior' is required for brma models.")
  }
  draws_df <- posterior::as_draws_df(model$fit$mcmc)
  mu <- as.numeric(draws_df[["mu"]])
  alpha <- (1 - prob) / 2
  c(estimate = stats::median(mu),
    lower    = unname(stats::quantile(mu, alpha)),
    upper    = unname(stats::quantile(mu, 1 - alpha)))
}

#' @export
extract_forest_draws.brma <- function(model, data, ...) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    rlang::abort("Package 'posterior' is required for brma models.")
  }
  if (!requireNamespace("RoBMA", quietly = TRUE)) {
    rlang::abort("Package 'RoBMA' is required.")
  }
  if (!"Author" %in% names(data)) {
    rlang::abort("Data must contain an 'Author' column.")
  }

  draws_df    <- posterior::as_draws_df(model$fit$mcmc)
  mu_draws    <- as.numeric(draws_df[["mu"]])
  n_draws     <- length(mu_draws)
  draw_ids    <- seq_len(n_draws)
  study_names <- as.character(data$Author)
  n_studies   <- length(study_names)

  tau_vec <- if ("tau" %in% names(draws_df)) {
    as.numeric(draws_df[["tau"]])
  } else {
    rep(NA_real_, n_draws)
  }

  blup_mat <- tryCatch({
    unclass(stats::predict(model, newdata = NULL, type = "estimate", quiet = TRUE))
  }, error = function(e) NULL)

  study_draws <- purrr::map_dfr(seq_len(n_studies), function(i) {
    theta_i <- if (!is.null(blup_mat) && i <= ncol(blup_mat)) {
      as.numeric(blup_mat[, i])
    } else {
      mu_draws
    }

    tibble::tibble(
      .draw       = draw_ids,
      Author      = study_names[i],
      theta_study = theta_i,
      eta         = theta_i - mu_draws,
      tau         = tau_vec
    )
  })

  pooled_draws <- tibble::tibble(
    .draw       = draw_ids,
    Author      = "Pooled Effect",
    theta_study = mu_draws,
    eta         = NA_real_,
    tau         = tau_vec
  )

  dplyr::bind_rows(study_draws, pooled_draws)
}

#' @export
refit_for_subgroup.brma <- function(model, data, ...) {
  if (!requireNamespace("RoBMA", quietly = TRUE)) {
    rlang::abort("Package 'RoBMA' is required.")
  }

  tryCatch({
    outcome_priors <- model$priors$outcome
    fit_ctrl       <- model$fit_control
    measure        <- attr(model$data, "measure")

    if (!all(c("yi", "vi") %in% names(data))) {
      rlang::abort("Data must contain 'yi' and 'vi' columns for brma refit.")
    }

    sg_data <- data.frame(
      yi   = data$yi,
      sei  = sqrt(data$vi),
      slab = as.character(data$Author)
    )

    refit_fn <- if (inherits(model, "RoBMA")) RoBMA::RoBMA else RoBMA::brma

    args <- list(
      yi      = sg_data$yi,
      sei     = sg_data$sei,
      slab    = sg_data$slab,
      measure = measure,
      sample  = fit_ctrl$sample  %||% 5000L,
      burnin  = fit_ctrl$burnin  %||% 2000L,
      adapt   = fit_ctrl$adapt   %||% 500L,
      chains  = fit_ctrl$chains  %||% 3L,
      thin    = fit_ctrl$thin    %||% 1L,
      seed    = fit_ctrl$seed,
      silent  = TRUE
    )

    if (!is.null(outcome_priors$mu)) {
      args$prior_effect <- outcome_priors$mu
    }
    if (!is.null(outcome_priors$tau)) {
      args$prior_heterogeneity <- outcome_priors$tau
    }

    if (inherits(model, "RoBMA") && !is.null(outcome_priors$bias)) {
      args$prior_bias <- outcome_priors$bias
    }

    do.call(refit_fn, args)
  }, error = function(e) {
    rlang::warn(
      c(
        paste0("brma subgroup refit failed: ", e$message),
        "i" = "Subgroup pooled effects will be approximated from the overall model."
      )
    )
    NULL
  })
}


# =============================================================================
# JAGS (R2jags / rjags) METHODS
# =============================================================================

#' @export
detect_model_info.rjags <- function(model, tau_par = "tau", ...) {
  sims <- model$BUGSoutput$sims.matrix
  tau_found <- any(grepl(
    paste0("^", tau_par, "$|^sigma$|^sd$|^tau$"),
    colnames(sims)
  ))
  list(backend = "jags",
       pooling = if (tau_found) "random_effect" else "common_effect")
}

#' @export
extract_pooled_summary.rjags <- function(model, prob = 0.95,
                                         theta_par = "mu", ...) {
  sims <- model$BUGSoutput$sims.matrix
  theta_col <- .find_parameter(as.data.frame(sims), theta_par,
                               c("theta", "d", "beta"))
  theta <- sims[, theta_col]
  alpha <- (1 - prob) / 2
  c(estimate = stats::median(theta),
    lower    = unname(stats::quantile(theta, alpha)),
    upper    = unname(stats::quantile(theta, 1 - alpha)))
}

#' @export
extract_forest_draws.rjags <- function(model, data,
                                       theta_par = "mu",
                                       tau_par = "sigma",
                                       study_par = NULL, ...) {
  if (!"Author" %in% names(data)) {
    rlang::abort("Data must contain an 'Author' column.")
  }
  posterior_df <- as.data.frame(model$BUGSoutput$sims.matrix)
  .extract_stan_draws(posterior_df, data, theta_par, tau_par, study_par)
}

#' @export
refit_for_subgroup.rjags <- function(model, data, original_data = NULL, ...) {
  if (!requireNamespace("R2jags", quietly = TRUE)) {
    rlang::abort("Package 'R2jags' is required to refit JAGS models.")
  }

  old_data <- model$model$data()
  new_data <- .rebuild_jags_data(old_data, data, original_data)
  if (is.null(new_data)) return(NULL)

  tryCatch({
    R2jags::jags(
      data               = new_data,
      parameters.to.save = model$parameters.to.save,
      model.file         = model$model.file,
      n.iter             = model$BUGSoutput$n.iter,
      n.burnin           = model$BUGSoutput$n.burnin,
      n.thin             = model$BUGSoutput$n.thin,
      n.chains           = model$BUGSoutput$n.chains
    )
  },
  error = function(e) {
    rlang::warn(
      c("JAGS refit failed for subgroup.",
        "i" = conditionMessage(e)),
      .frequency = "once",
      .frequency_id = "refit_jags_error"
    )
    NULL
  })
}


#' Rebuild JAGS data list for a subgroup
#'
#' Subsets array/vector elements in the original JAGS data list to match the
#' rows present in `data`. Scalars and matrices without a study dimension are
#' kept as-is.
#'
#' @noRd
.rebuild_jags_data <- function(old_data, data, original_data) {
  n_old <- old_data[["S"]] %||% old_data[["N"]] %||% old_data[["n.studies"]]
  if (is.null(n_old)) {
    rlang::warn(
      c("Cannot detect study count variable in JAGS data.",
        "i" = "Expected 'S', 'N', or 'n.studies' in the data list."),
      .frequency = "once",
      .frequency_id = "refit_jags_no_S"
    )
    return(NULL)
  }

  rows <- match(data$Author, original_data$Author)
  rows <- rows[!is.na(rows)]
  if (length(rows) == 0) return(NULL)

  n_new <- length(rows)
  new_data <- old_data

  for (nm in names(old_data)) {
    val <- old_data[[nm]]
    if (is.vector(val) && length(val) == n_old) {
      new_data[[nm]] <- val[rows]
    } else if (is.matrix(val) && nrow(val) == n_old) {
      new_data[[nm]] <- val[rows, , drop = FALSE]
    }
  }

  study_count_name <- if ("S" %in% names(old_data)) "S"
    else if ("N" %in% names(old_data)) "N"
    else "n.studies"
  new_data[[study_count_name]] <- n_new

  new_data
}


# =============================================================================
# INTERNAL HELPERS
# =============================================================================

#' Find a parameter column in a posterior data frame
#' @noRd
.find_parameter <- function(posterior_df, primary, fallbacks = character(0)) {
  col <- grep(paste0("^", primary, "$|^", primary, "\\["),
              names(posterior_df), value = TRUE)
  if (length(col) > 0) return(col[1])

  for (fb in fallbacks) {
    col <- grep(paste0("^", fb, "$|^", fb, "\\["),
                names(posterior_df), value = TRUE)
    if (length(col) > 0) return(col[1])
  }

  rlang::abort(
    c(
      paste0("Could not find parameter '", primary, "' in posterior samples."),
      "i" = paste0("Available: ",
                    paste(head(names(posterior_df), 20), collapse = ", "))
    )
  )
}

#' Shared Stan draw extraction for stanfit and CmdStanMCMC
#' @noRd
.extract_stan_draws <- function(posterior_df, data,
                                theta_par, tau_par, study_par) {
  n_studies   <- nrow(data)
  study_names <- as.character(data$Author)
  n_draws     <- nrow(posterior_df)
  draw_ids    <- seq_len(n_draws)

  theta_col <- .find_parameter(posterior_df, theta_par, c("mu", "beta", "d"))
  theta     <- posterior_df[[theta_col]]

  tau_col <- grep(paste0("^", tau_par, "$|^", tau_par, "\\[1\\]$"),
                  names(posterior_df), value = TRUE)
  if (length(tau_col) == 0) {
    tau_col <- grep("^sigma$|^sd$", names(posterior_df), value = TRUE)
  }
  tau_vec <- if (length(tau_col) > 0) posterior_df[[tau_col[1]]] else rep(NA_real_, n_draws)

  # Auto-detect study-level parameterisation
  effect_type <- "none"
  if (is.null(study_par)) {
    candidates <- c(eta = "deviation", u = "deviation",
                    mu = "absolute", d = "absolute")
    for (nm in names(candidates)) {
      if (any(grepl(paste0("^", nm, "\\["), names(posterior_df)))) {
        study_par   <- nm
        effect_type <- candidates[[nm]]
        break
      }
    }
  } else {
    effect_type <- if (study_par %in% c("eta", "u", "eps", "epsilon")) {
      "deviation"
    } else {
      "absolute"
    }
  }

  study_draws <- purrr::map_dfr(seq_len(n_studies), function(i) {
    if (effect_type == "none" || is.null(study_par)) {
      deviation    <- rep(0, n_draws)
      study_effect <- theta
    } else {
      effect_col <- paste0(study_par, "[", i, "]")
      if (effect_col %in% names(posterior_df)) {
        if (effect_type == "deviation") {
          deviation    <- posterior_df[[effect_col]]
          study_effect <- theta + deviation
        } else {
          study_effect <- posterior_df[[effect_col]]
          deviation    <- study_effect - theta
        }
      } else {
        deviation    <- rep(0, n_draws)
        study_effect <- theta
      }
    }

    tibble::tibble(
      .draw       = draw_ids,
      Author      = study_names[i],
      theta_study = study_effect,
      eta         = deviation,
      tau         = tau_vec
    )
  })

  pooled_draws <- tibble::tibble(
    .draw       = draw_ids,
    Author      = "Pooled Effect",
    theta_study = theta,
    eta         = NA_real_,
    tau         = tau_vec
  )

  dplyr::bind_rows(study_draws, pooled_draws)
}
