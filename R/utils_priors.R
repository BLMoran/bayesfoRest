
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
#'   \item normal() to 𝒩()
#'   \item cauchy() to 𝒞()
#'   \item student_t() to 𝓉()
#'   \item exponential() to ℰ()
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
  
  prior_string <- trimws(prior_string)
  
  # Normal
  if (grepl("^normal\\s*\\(", prior_string)) {
    params <- sub("^normal\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    params <- gsub("\\s+", "", params)
    # Add a thin space after closing paren - if repeated, won't be visible
    return(paste0(script_N, "(", params, ")\u2009"))
  }
  
  # Cauchy (ignore truncation in display)
  if (grepl("^cauchy\\s*\\(", prior_string)) {
    base <- sub("\\s*\\[.*$", "", prior_string)
    params <- sub("^cauchy\\s*\\(([^)]+)\\).*$", "\\1", base)
    params <- gsub("\\s+", "", params)
    return(paste0(script_C, "(", params, ")\u2009"))
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

#' Extract Mu and Tau Priors
#'
#' @description
#' Extracts and formats the prior specifications for both the mean (mu) and
#' heterogeneity (tau) parameters from a list of priors.
#'
#' @param priors Named list of prior specifications
#' @param model A fitted brmsfit object (for extracting default tau prior)
#'
#' @return A tibble with columns: prior, prior_label, mu_prior_unicode, tau_prior_unicode
#'
#' @keywords internal
#' @noRd
extract_mu_tau_priors <- function(priors, model) {
  
  model_tau <- extract_tau_from_model(model)
  
  purrr::imap_dfr(priors, function(prior_obj, prior_name) {
    
    mu_prior <- prior_obj |>
      dplyr::filter(class == "Intercept") |>
      dplyr::pull(prior)
    
    tau_prior <- prior_obj |>
      dplyr::filter(class == "sd") |>
      dplyr::pull(prior)
    
    tibble::tibble(
      prior = prior_name,
      prior_label = dplyr::recode(
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
#' Extracts and formats prior specifications from a RoBMA_brms fit object
#' for display in tables.
#'
#' @param robma_fit A RoBMA_brms fit object
#'
#' @return List with mu_prior_unicode and tau_prior_unicode
#'
#' @keywords internal
#' @noRd
extract_priors_from_robma_fit <- function(robma_fit) {
  
  # ---- μ priors ----
  mu_strings <- robma_fit$prior_specs$mu$mu_prior_string
  
  mu_prior_unicode <- paste(
    purrr::map_chr(
      mu_strings,
      ~ paste0(prior_to_unicode(.x), "\u2009")
    ),
    collapse = " / "
  )
  
  # ---- τ priors ----
  tau_strings <- robma_fit$prior_specs$tau$tau_prior_string
  
  tau_prior_unicode <- paste(
    purrr::map_chr(
      tau_strings,
      function(x) {
        if (x == "fixed") {
          paste0("τ = 0", "\u2009")
        } else {
          paste0(prior_to_unicode(x), "\u2009")
        }
      }
    ),
    collapse = " / "
  )
  
  list(
    mu_prior_unicode  = mu_prior_unicode,
    tau_prior_unicode = tau_prior_unicode
  )
}
