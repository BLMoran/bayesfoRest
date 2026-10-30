
#' Extract Tau Prior from Model
#'
#' @description
#' Extracts the prior specification for the heterogeneity parameter (tau)
#' from a brms model.
#'
#' @param model A fitted brmsfit object
#'
#' @return Character string of the tau prior, or NA if not found
#'
#' @keywords internal
#' @noRd
extract_tau_from_model <- function(model) {

  if (!inherits(model, "brmsfit")) return(NA_character_)

  tau_prior <- model$prior |>
    dplyr::filter(class == "sd") |>
    dplyr::pull(prior)

  if (length(tau_prior) == 0) {
    NA_character_
  } else {
    tau_prior[[1]]
  }
}

#' Convert Prior String to Unicode Format
#'
#' @description
#' Converts prior distribution strings to Unicode mathematical notation
#' for pretty printing in tables.
#'
#' @param prior_string Character string representing a prior distribution
#'
#' @return Unicode-formatted prior string
#'
#' @details
#' Supports conversion of:
#' \itemize{
#'   \item normal() to \eqn{\mathcal{N}()}
#'   \item cauchy() to \eqn{\mathcal{C}()}
#'   \item student_t() to \eqn{\mathcal{t}()}
#'   \item exponential() to \eqn{\mathcal{E}()}
#' }
#'
#' @keywords internal
#' @noRd
prior_to_unicode <- function(prior_string) {
  if (is.null(prior_string) || is.na(prior_string)) {
    return(NA_character_)
  }

  # Unicode math-script letters (STIX Two Math)
  script_N <- "\U0001D4A9"  # 𝒩
  script_C <- "\U0001D49E"  # 𝒞
  script_t <- "\U0001D4C9"  # 𝓉
  script_E <- "\U2130"      # ℰ
  script_I <- "\U2110"      # ℐ
  script_G <- "\U0001D4A2"  # 𝒢
  delta    <- "\u03B4"      # δ

  prior_string <- trimws(prior_string)

  # Point/spike prior → δ(0) (Dirac delta)
  if (grepl("^point\\s*\\(", prior_string)) {
    loc <- sub("^point\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    return(paste0(delta, "(", loc, ")\u2009"))
  }

  # Inverse gamma → ℐ𝒢(shape, scale)
  if (grepl("^inv_gamma\\s*\\(", prior_string)) {
    params <- sub("^inv_gamma\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    params <- gsub("\\s+", "", params)
    return(paste0(script_I, script_G, "(", params, ")\u2009"))
  }

  # Normal
  if (grepl("^normal\\s*\\(", prior_string)) {
    params <- sub("^normal\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    params <- gsub("\\s+", "", params)
    return(paste0(script_N, "(", params, ")\u2009"))
  }

  # Half-cauchy
  if (grepl("^half_cauchy\\s*\\(", prior_string)) {
    params <- sub("^half_cauchy\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    params <- gsub("\\s+", "", params)
    return(paste0(script_C, "\u207A(", params, ")\u2009"))
  }

  # Cauchy (ignore truncation in display)
  if (grepl("^cauchy\\s*\\(", prior_string)) {
    base <- sub("\\s*\\[.*$", "", prior_string)
    params <- sub("^cauchy\\s*\\(([^)]+)\\).*$", "\\1", base)
    params <- gsub("\\s+", "", params)
    return(paste0(script_C, "(", params, ")\u2009"))
  }

  # Half-normal
  if (grepl("^half_normal\\s*\\(", prior_string)) {
    params <- sub("^half_normal\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    params <- gsub("\\s+", "", params)
    return(paste0(script_N, "\u207A(", params, ")\u2009"))
  }

  # Student-t
  if (grepl("^student_t\\s*\\(", prior_string)) {
    clean <- sub("\\[.*$", "", prior_string)
    params <- sub("^student_t\\s*\\(([^)]+)\\).*$", "\\1", clean)
    parts <- strsplit(params, ",")[[1]]
    parts <- trimws(parts)

    if (length(parts) != 3) {
      return(prior_string)
    }

    df  <- gsub("\\s+", "", parts[1])
    loc <- gsub("\\s+", "", parts[2])
    scl <- gsub("\\s+", "", parts[3])

    # Subscript df if single digit
    if (nchar(df) == 1 && grepl("^[0-9]$", df)) {
      df <- c(
        "0" = "\u2080", "1" = "\u2081", "2" = "\u2082", "3" = "\u2083",
        "4" = "\u2084", "5" = "\u2085", "6" = "\u2086", "7" = "\u2087",
        "8" = "\u2088", "9" = "\u2089"
      )[df]
    }

    return(paste0(script_t, df, "(", loc, ", ", scl, ")\u2009"))
  }

  # Exponential
  if (grepl("^exponential\\s*\\(", prior_string)) {
    params <- sub("^exponential\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    params <- gsub("\\s+", "", params)
    return(paste0(script_E, "(", params, ")\u2009"))
  }

  prior_string
}

#' Extract Mu and Tau Priors from Sensitivity Prior List
#'
#' @description
#' Extracts and formats the prior specifications for both the mean (mu) and
#' heterogeneity (tau) parameters from a named list of priors. Supports
#' \code{bf_prior} objects (the standard format) as well as legacy brms
#' prior objects.
#'
#' @param priors Named list of prior specifications. Each element should be a
#'   list with \code{name}, \code{mu_prior}, and \code{tau_prior}.
#' @param model A fitted model object (used for extracting default tau prior
#'   when using legacy brms format).
#'
#' @return A tibble with columns: prior, prior_label, mu_prior_unicode,
#'   tau_prior_unicode
#'
#' @keywords internal
#' @noRd
extract_mu_tau_priors <- function(priors, model) {

  purrr::imap_dfr(priors, function(prior_spec, prior_name) {

    # --- New bf_prior format ---
    if (inherits(prior_spec$mu_prior, "bf_prior") ||
        inherits(prior_spec$tau_prior, "bf_prior")) {

      mu_unicode  <- format_bf_prior(prior_spec$mu_prior)
      tau_unicode <- format_bf_prior(prior_spec$tau_prior)

      return(tibble::tibble(
        prior             = prior_name,
        prior_label       = prior_spec$name %||% prior_name,
        mu_prior_unicode  = mu_unicode,
        tau_prior_unicode = tau_unicode
      ))
    }

    # --- Legacy brms prior format ---
    if (inherits(prior_spec, "brmsprior")) {
      model_tau <- extract_tau_from_model(model)

      mu_prior <- prior_spec |>
        dplyr::filter(class == "Intercept") |>
        dplyr::pull(prior)

      tau_prior <- prior_spec |>
        dplyr::filter(class == "sd") |>
        dplyr::pull(prior)

      return(tibble::tibble(
        prior = prior_name,
        prior_label = prior_spec$name %||%
          dplyr::recode(
            prior_name,
            vague = "Vague",
            weakreg = "Weakly Regularising",
            informative = "Informative",
            .default = stringr::str_to_sentence(prior_name)
          ),
        mu_prior_unicode =
          if (length(mu_prior) > 0) prior_to_unicode(mu_prior[[1]]) else NA_character_,
        tau_prior_unicode =
          if (length(tau_prior) > 0) {
            prior_to_unicode(tau_prior[[1]])
          } else if (!is.na(model_tau)) {
            prior_to_unicode(model_tau)
          } else {
            NA_character_
          }
      ))
    }

    # --- Fallback: list with name only (no bf_prior objects) ---
    tibble::tibble(
      prior             = prior_name,
      prior_label       = prior_spec$name %||% prior_name,
      mu_prior_unicode  = NA_character_,
      tau_prior_unicode = NA_character_
    )
  })
}

#' Apply Math Font to GT Table
#'
#' @description
#' Applies a mathematical font (like STIX Two Math) to specific columns
#' in a gt table for proper rendering of mathematical symbols.
#'
#' @param gt_tbl A gt table object
#' @param columns Column names to apply the math font to
#' @param math_font Character string. Name of the math font. Default is "STIX Two Math"
#'
#' @return Modified gt table object
#'
#' @keywords internal
#' @noRd
apply_math_font <- function(
    gt_tbl,
    columns,
    math_font = "STIX Two Math"
) {
  gt_tbl |>
    gt::tab_style(
      style = gt::cell_text(font = math_font),
      locations = gt::cells_body(columns = columns))
}

#' Extract Priors from RoBMA Fit
#'
#' @description
#' Extracts and formats prior specifications from a fitted RoBMA object
#' for display in sensitivity analysis tables. Supports both native RoBMA
#' objects and the legacy RoBMA_brms wrapper format.
#'
#' @param robma_fit A fitted RoBMA object or RoBMA_brms wrapper
#'
#' @return List with mu_prior_unicode and tau_prior_unicode
#'
#' @keywords internal
#' @noRd
extract_priors_from_robma_fit <- function(robma_fit) {

  # --- Legacy RoBMA_brms wrapper format ---
  if (!is.null(robma_fit$prior_specs)) {
    mu_strings <- robma_fit$prior_specs$mu$mu_prior_string
    mu_prior_unicode <- paste(
      purrr::map_chr(mu_strings, ~ paste0(prior_to_unicode(.x), "\u2009")),
      collapse = " / "
    )

    tau_strings <- robma_fit$prior_specs$tau$tau_prior_string
    tau_prior_unicode <- paste(
      purrr::map_chr(tau_strings, function(x) {
        if (x == "fixed") paste0("\u03C4 = 0", "\u2009")
        else paste0(prior_to_unicode(x), "\u2009")
      }),
      collapse = " / "
    )

    return(list(
      mu_prior_unicode  = mu_prior_unicode,
      tau_prior_unicode = tau_prior_unicode
    ))
  }

  # --- Native RoBMA object ---
  # Extract unique mu and tau priors from the ensemble's models
  models <- robma_fit[["models"]]
  if (is.null(models)) {
    return(list(mu_prior_unicode = NA_character_, tau_prior_unicode = NA_character_))
  }

  mu_strings  <- character()
  tau_strings <- character()

  for (m in models) {
    prior_names <- names(m[["priors"]])
    for (nm in prior_names) {
      p <- m[["priors"]][[nm]]
      prior_str <- .robma_prior_to_string(p)

      if (nm == "mu" && !prior_str %in% mu_strings) {
        mu_strings <- c(mu_strings, prior_str)
      } else if (nm == "tau" && !prior_str %in% tau_strings) {
        tau_strings <- c(tau_strings, prior_str)
      }
    }
  }

  mu_prior_unicode <- if (length(mu_strings) > 0) {
    paste(purrr::map_chr(mu_strings, ~ prior_to_unicode(.x)), collapse = " / ")
  } else {
    NA_character_
  }

  tau_prior_unicode <- if (length(tau_strings) > 0) {
    paste(purrr::map_chr(tau_strings, function(x) {
      if (x == "point(0)") "\u03C4 = 0" else prior_to_unicode(x)
    }), collapse = " / ")
  } else {
    NA_character_
  }

  list(
    mu_prior_unicode  = mu_prior_unicode,
    tau_prior_unicode = tau_prior_unicode
  )
}


#' Convert RoBMA Prior Object to String
#'
#' @description
#' Converts a BayesTools prior object from a RoBMA model to a human-readable
#' distribution string (e.g. "normal(0,1)").
#'
#' @param p A prior object from a RoBMA model.
#'
#' @return Character string representing the prior distribution.
#'
#' @keywords internal
#' @noRd
.robma_prior_to_string <- function(p) {

  # Point/spike prior
  if (inherits(p, "prior.point") || inherits(p, "prior.spike")) {
    loc <- p[["parameters"]][["location"]] %||% 0
    return(paste0("point(", loc, ")"))
  }

  # No prior (none)
  if (inherits(p, "prior.none")) {
    return("none")
  }

  # Distribution priors — extract distribution name and parameters
  dist   <- p[["distribution"]]
  params <- p[["parameters"]]
  trunc  <- p[["truncation"]]

  if (is.null(dist) || is.null(params)) {
    return(tryCatch(
      paste0(capture.output(print(p)), collapse = ""),
      error = function(e) "unknown"
    ))
  }

  # Detect half-distribution (truncated at 0 from below)
  is_half <- !is.null(trunc) &&
    is.finite(trunc$lower %||% -Inf) &&
    (trunc$lower %||% -Inf) == 0 &&
    !is.finite(trunc$upper %||% Inf)

  # RoBMA stores cauchy as t(df=1); map back to recognizable names
  if (dist == "t" && !is.null(params$df) && params$df == 1) {
    loc <- params$location %||% 0
    scl <- params$scale %||% 1
    prefix <- if (is_half) "half_cauchy" else "cauchy"
    return(paste0(prefix, "(", loc, ",", scl, ")"))
  }

  if (dist == "t") {
    loc <- params$location %||% 0
    scl <- params$scale %||% 1
    df  <- params$df %||% 3
    prefix <- if (is_half) "half_student_t" else "student_t"
    return(paste0(prefix, "(", df, ",", loc, ",", scl, ")"))
  }

  if (dist == "normal") {
    mn <- params$mean %||% 0
    sd <- params$sd %||% 1
    prefix <- if (is_half) "half_normal" else "normal"
    return(paste0(prefix, "(", mn, ",", sd, ")"))
  }

  if (dist == "invgamma") {
    sh <- params$shape %||% 1
    sc <- params$scale %||% 1
    return(paste0("inv_gamma(", sh, ",", sc, ")"))
  }

  # Fallback: generic
  param_str <- paste(unlist(params), collapse = ",")
  prefix <- if (is_half) paste0("half_", dist) else dist
  paste0(prefix, "(", param_str, ")")
}
