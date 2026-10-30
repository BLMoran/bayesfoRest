#' Prior Distribution Constructors for Sensitivity Analysis
#'
#' @description
#' Backend-agnostic prior specification functions for use with
#' \code{\link{sensitivity_plot}}. These create lightweight prior objects
#' that are automatically translated to the appropriate format for each
#' model backend (brms, bayesmeta, MetaStan, rstan, cmdstanr).
#'
#' @param mean,location Prior mean or location parameter.
#' @param sd,scale Prior standard deviation or scale parameter.
#' @param df Degrees of freedom (for \code{half_t} only).
#' @param rate Rate parameter (for \code{bf_exponential} only).
#'
#' @return An object of class \code{"bf_prior"} containing the distribution
#'   family and parameters.
#'
#' @details
#' The following distributions are available:
#' \describe{
#'   \item{\code{bf_normal(mean, sd)}}{Normal distribution. Typically used for
#'     the overall effect (\eqn{\mu}) prior.}
#'   \item{\code{bf_half_normal(mean, sd)}}{Half-normal distribution (truncated
#'     at 0). Typically used for the heterogeneity (\eqn{\tau}) prior.}
#'   \item{\code{bf_half_cauchy(location, scale)}}{Half-Cauchy distribution.
#'     A common default for \eqn{\tau} priors.}
#'   \item{\code{bf_half_t(df, location, scale)}}{Half-Student-t distribution.
#'     Flexible heavy-tailed prior for \eqn{\tau}.}
#'   \item{\code{bf_exponential(rate)}}{Exponential distribution for
#'     \eqn{\tau}.}
#' }
#'
#' @examples
#' \dontrun{
#' priors <- list(
#'   vague = list(
#'     name      = "Vague",
#'     mu_prior  = bf_normal(0, 10),
#'     tau_prior = bf_half_cauchy(0, 1)
#'   ),
#'   weakreg = list(
#'     name      = "Weakly Regularising",
#'     mu_prior  = bf_normal(0, 1),
#'     tau_prior = bf_half_cauchy(0, 0.5)
#'   ),
#'   informative = list(
#'     name      = "Informative",
#'     mu_prior  = bf_normal(0, 0.5),
#'     tau_prior = bf_half_cauchy(0, 0.3)
#'   )
#' )
#'
#' sensitivity_plot(model, data, priors = priors, measure = "OR")
#' }
#'
#' @name bf_priors
NULL

#' @rdname bf_priors
#' @export
bf_normal <- function(mean = 0, sd = 1) {
  structure(
    list(family = "normal", mean = mean, sd = sd),
    class = "bf_prior"
  )
}

#' @rdname bf_priors
#' @export
bf_half_normal <- function(mean = 0, sd = 1) {
  structure(
    list(family = "half_normal", mean = mean, sd = sd),
    class = "bf_prior"
  )
}

#' @rdname bf_priors
#' @export
bf_half_cauchy <- function(location = 0, scale = 1) {
  structure(
    list(family = "half_cauchy", location = location, scale = scale),
    class = "bf_prior"
  )
}

#' @rdname bf_priors
#' @export
bf_half_t <- function(df = 3, location = 0, scale = 1) {
  structure(
    list(family = "half_t", df = df, location = location, scale = scale),
    class = "bf_prior"
  )
}

#' @rdname bf_priors
#' @export
bf_exponential <- function(rate = 1) {
  structure(
    list(family = "exponential", rate = rate),
    class = "bf_prior"
  )
}

#' @export
print.bf_prior <- function(x, ...) {
  cat(format_bf_prior(x), "\n")
}


# =========================================================================
# Translation to backend-specific formats
# =========================================================================

#' Translate bf_prior to brms Prior String
#'
#' @param prior A \code{bf_prior} object.
#' @param class Character. brms prior class ("Intercept" for mu, "sd" for tau).
#'
#' @return A \code{brmsprior} object.
#'
#' @keywords internal
#' @noRd
bf_prior_to_brms <- function(prior, class = "Intercept") {
  if (!requireNamespace("brms", quietly = TRUE)) {
    rlang::abort("Package 'brms' is required.")
  }

  str <- switch(
    prior$family,
    normal      = paste0("normal(", prior$mean, ", ", prior$sd, ")"),
    half_normal = paste0("normal(", prior$mean, ", ", prior$sd, ")"),
    half_cauchy = paste0("cauchy(", prior$location, ", ", prior$scale, ")"),
    half_t      = paste0("student_t(", prior$df, ", ", prior$location, ", ", prior$scale, ")"),
    exponential = paste0("exponential(", prior$rate, ")"),
    rlang::abort(paste0("Unknown prior family: ", prior$family))
  )

  # brms handles truncation via class (sd class is automatically lower-bounded at 0)
  brms::prior_string(str, class = class)
}


#' Translate bf_prior to bayesmeta Prior Function
#'
#' @param prior A \code{bf_prior} object.
#' @param type Character. "mu" or "tau".
#'
#' @return A function suitable for bayesmeta's mu.prior or tau.prior argument.
#'
#' @keywords internal
#' @noRd
bf_prior_to_bayesmeta <- function(prior, type = c("mu", "tau")) {
  type <- match.arg(type)

  if (!requireNamespace("bayesmeta", quietly = TRUE)) {
    rlang::abort("Package 'bayesmeta' is required.")
  }

  switch(
    prior$family,

    normal = {
      m <- prior$mean; s <- prior$sd
      if (type == "mu") {
        # bayesmeta mu.prior: c(mean, sd)
        c(m, s)
      } else {
        # tau prior as density function (normal truncated at 0)
        function(t) dnorm(t, mean = m, sd = s)
      }
    },

    half_normal = {
      m <- prior$mean; s <- prior$sd
      if (type == "mu") {
        c(m, s)
      } else {
        function(t) bayesmeta::dhalfnormal(t, scale = s)
      }
    },

    half_cauchy = {
      loc <- prior$location; sc <- prior$scale
      if (type == "mu") {
        # Use normal approximation for mu (Cauchy not standard for mu)
        c(loc, sc)
      } else {
        function(t) bayesmeta::dhalfcauchy(t, scale = sc)
      }
    },

    half_t = {
      nu <- prior$df; loc <- prior$location; sc <- prior$scale
      if (type == "mu") {
        c(loc, sc)
      } else {
        function(t) bayesmeta::dhalft(t, df = nu, scale = sc)
      }
    },

    exponential = {
      r <- prior$rate
      if (type == "mu") {
        rlang::warn("Exponential prior is not appropriate for mu; using normal(0, 1/rate) instead.")
        c(0, 1 / r)
      } else {
        function(t) dexp(t, rate = r)
      }
    },

    rlang::abort(paste0("Unknown prior family: ", prior$family))
  )
}


#' Build Backend-Specific Prior Set from bf_prior Specifications
#'
#' @description
#' Translates a list containing \code{mu_prior} and \code{tau_prior}
#' (both \code{bf_prior} objects) into the format expected by a specific
#' model backend.
#'
#' @param prior_spec A list with elements \code{mu_prior} and \code{tau_prior},
#'   each a \code{bf_prior} object.
#' @param model A fitted model object (used to determine backend).
#'
#' @return Backend-specific prior object suitable for \code{update_model_prior()}.
#'
#' @keywords internal
#' @noRd
translate_priors_for_backend <- function(prior_spec, model) {

  mu_prior  <- prior_spec$mu_prior
  tau_prior <- prior_spec$tau_prior

  if (inherits(model, "brmsfit")) {
    # Combine into a single brmsprior data frame
    mu_brms  <- bf_prior_to_brms(mu_prior, class = "Intercept")
    tau_brms <- if (!is.null(tau_prior)) {
      bf_prior_to_brms(tau_prior, class = "sd")
    } else {
      NULL
    }
    if (!is.null(tau_brms)) {
      rbind(mu_brms, tau_brms)
    } else {
      mu_brms
    }

  } else if (inherits(model, "bayesmeta")) {
    list(
      mu_prior  = bf_prior_to_bayesmeta(mu_prior, "mu"),
      tau_prior = if (!is.null(tau_prior)) {
        bf_prior_to_bayesmeta(tau_prior, "tau")
      } else {
        NULL
      }
    )

  } else if (inherits(model, "meta_stan")) {
    # MetaStan uses c(location, scale) numeric vectors
    # tau_prior_dist tells MetaStan which distribution ("half-normal", "half-cauchy", etc.)
    mu_metastan <- switch(
      mu_prior$family,
      normal      = c(mu_prior$mean, mu_prior$sd),
      half_normal = c(mu_prior$mean, mu_prior$sd),
      half_cauchy = c(mu_prior$location, mu_prior$scale),
      half_t      = c(mu_prior$location, mu_prior$scale),
      c(0, 10)
    )

    # MetaStan tau_prior is a SCALAR (the scale/rate parameter), not a vector
    tau_metastan <- NULL
    tau_dist     <- NULL
    if (!is.null(tau_prior)) {
      tau_metastan <- switch(
        tau_prior$family,
        half_normal = tau_prior$sd,
        half_cauchy = tau_prior$scale,
        half_t      = tau_prior$scale,
        exponential = tau_prior$rate,
        1
      )
      tau_dist <- switch(
        tau_prior$family,
        half_normal = "half-normal",
        half_cauchy = "half-cauchy",
        half_t      = "half-t",
        exponential = "exponential",
        "half-normal"
      )
    }

    list(
      mu_prior       = mu_metastan,
      tau_prior      = tau_metastan,
      tau_prior_dist = tau_dist
    )

  } else {
    # stanfit / CmdStanMCMC / stanreg — cannot easily update priors
    # Return the raw bf_prior objects; update_model_prior will handle
    list(
      mu_prior  = mu_prior,
      tau_prior = tau_prior
    )
  }
}


# =========================================================================
# Formatting for display
# =========================================================================

#' Format bf_prior as Unicode String
#'
#' @param prior A \code{bf_prior} object.
#'
#' @return Character string with unicode math symbols.
#'
#' @keywords internal
#' @noRd
format_bf_prior <- function(prior) {

  if (is.null(prior)) return(NA_character_)
  if (!inherits(prior, "bf_prior")) return(NA_character_)

  # Unicode math-script letters
  script_N <- "\U0001D4A9"  # 𝒩
  script_C <- "\U0001D49E"  # 𝒞
  script_t <- "\U0001D4C9"  # 𝓉
  script_E <- "\U2130"      # ℰ
  superscript_plus <- "\u207A"  # ⁺

  # Trailing thin space (\u2009) prevents STIX Two Math font from

  # rendering a doubled closing parenthesis in gt tables.
  thin <- "\u2009"

  switch(
    prior$family,

    normal = paste0(
      script_N, "(", prior$mean, ",", prior$sd, ")", thin
    ),

    half_normal = paste0(
      script_N, superscript_plus, "(", prior$mean, ",", prior$sd, ")", thin
    ),

    half_cauchy = paste0(
      script_C, superscript_plus, "(", prior$location, ",", prior$scale, ")", thin
    ),

    half_t = {
      df <- prior$df
      # Subscript digits
      sub_digits <- c(
        "0" = "\u2080", "1" = "\u2081", "2" = "\u2082", "3" = "\u2083",
        "4" = "\u2084", "5" = "\u2085", "6" = "\u2086", "7" = "\u2087",
        "8" = "\u2088", "9" = "\u2089"
      )
      df_str <- if (nchar(as.character(df)) == 1 && grepl("^[0-9]$", df)) {
        sub_digits[as.character(df)]
      } else {
        as.character(df)
      }
      paste0(
        script_t, df_str, superscript_plus,
        "(", prior$location, ",", prior$scale, ")", thin
      )
    },

    exponential = paste0(
      script_E, "(", prior$rate, ")", thin
    ),

    paste0("Unknown(", prior$family, ")")
  )
}
