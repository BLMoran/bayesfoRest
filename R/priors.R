#' Prior Distribution Constructors for Sensitivity Analysis
#'
#' @description
#' Backend-agnostic prior specification functions for use with
#' \code{\link{sensitivity_plot}}, \code{\link{ecdf_model_plot}},
#' \code{\link{ecdf_prior_plot}} and \code{\link{overall_plot}}. These create
#' lightweight prior objects that are translated to the format each model
#' backend expects. A prior a backend cannot fit exactly is an error, never a
#' silent substitution.
#'
#' @param mean,location Prior mean or location parameter.
#' @param sd,scale Prior standard deviation or scale parameter.
#' @param df Degrees of freedom.
#' @param rate Rate parameter (for \code{bf_exponential} only).
#' @param meanlog,sdlog Mean and standard deviation of \eqn{\log(\tau)}
#'   (for \code{bf_lognormal} only).
#' @param lower,upper Bounds of the uniform prior. For \eqn{\tau},
#'   \code{lower} must be 0.
#'
#' @return An object of class \code{"bf_prior"} containing the distribution
#'   family and parameters.
#'
#' @details
#' Priors for the pooled effect (\eqn{\mu}) are on the model scale (the log
#' scale for OR, RR, HR and IRR). Priors for the heterogeneity (\eqn{\tau})
#' are on the standard deviation of the study effects.
#'
#' \strong{Priors for \eqn{\mu}}
#' \describe{
#'   \item{\code{bf_normal(mean, sd)}}{Normal distribution.}
#'   \item{\code{bf_student_t(df, location, scale)}}{Student-t distribution:
#'     a heavier-tailed, weakly informative alternative to the normal.
#'     \code{df = 1} is a Cauchy distribution.}
#' }
#'
#' \strong{Priors for \eqn{\tau}}
#' \describe{
#'   \item{\code{bf_half_normal(mean, sd)}}{Half-normal distribution
#'     (truncated at 0).}
#'   \item{\code{bf_half_cauchy(location, scale)}}{Half-Cauchy distribution.}
#'   \item{\code{bf_half_t(df, location, scale)}}{Half-Student-t
#'     distribution.}
#'   \item{\code{bf_exponential(rate)}}{Exponential distribution.}
#'   \item{\code{bf_lognormal(meanlog, sdlog)}}{Log-normal distribution.
#'     Use \code{\link{bf_turner}} for the empirical priors of Turner et al.}
#'   \item{\code{bf_uniform(0, upper)}}{Uniform distribution from 0 to
#'     \code{upper}.}
#'   \item{\code{\link{bf_rhodes}()}}{Empirical log-t priors of Rhodes et al.}
#' }
#'
#' \strong{Backend support}
#'
#' Sensitivity refits stop with an error if the model's backend cannot fit a
#' prior exactly:
#' \tabular{lccc}{
#'   Prior \tab brms \tab bayesmeta \tab MetaStan \cr
#'   \code{bf_normal} (\eqn{\mu}) \tab yes \tab yes \tab yes \cr
#'   \code{bf_student_t} (\eqn{\mu}) \tab yes \tab no \tab no \cr
#'   \code{bf_half_normal} \tab yes \tab yes \tab yes (mean 0) \cr
#'   \code{bf_half_cauchy} \tab yes \tab yes \tab yes (location 0) \cr
#'   \code{bf_half_t} \tab yes \tab yes \tab no \cr
#'   \code{bf_exponential} \tab yes \tab yes \tab no \cr
#'   \code{bf_lognormal}, \code{bf_turner} \tab yes \tab yes \tab no \cr
#'   \code{bf_uniform} \tab yes \tab yes \tab yes \cr
#'   \code{bf_rhodes} \tab yes \tab yes \tab no \cr
#' }
#' Priors cannot be changed for rstan, cmdstanr, rstanarm, bayesma, RoBMA or
#' JAGS models.
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
#'     mu_prior  = bf_student_t(3, 0, 1),
#'     tau_prior = bf_half_normal(0, 0.5)
#'   ),
#'   empirical = list(
#'     name      = "Empirical (Turner)",
#'     mu_prior  = bf_normal(0, 1),
#'     tau_prior = bf_turner("all-cause mortality",
#'                           "pharmacological", "placebo / control")
#'   )
#' )
#'
#' sensitivity_plot(model, data, priors = priors, measure = "OR")
#' }
#'
#' @name bf_priors
#' @keywords internal
NULL

#' @rdname bf_priors
#' @export
bf_normal <- function(mean = 0, sd = 1) {
  new_bf_prior("normal", mean = mean, sd = sd)
}

#' @rdname bf_priors
#' @export
bf_student_t <- function(df = 3, location = 0, scale = 1) {
  new_bf_prior("student_t", df = df, location = location, scale = scale)
}

#' @rdname bf_priors
#' @export
bf_half_normal <- function(mean = 0, sd = 1) {
  new_bf_prior("half_normal", mean = mean, sd = sd)
}

#' @rdname bf_priors
#' @export
bf_half_cauchy <- function(location = 0, scale = 1) {
  new_bf_prior("half_cauchy", location = location, scale = scale)
}

#' @rdname bf_priors
#' @export
bf_half_t <- function(df = 3, location = 0, scale = 1) {
  new_bf_prior("half_t", df = df, location = location, scale = scale)
}

#' @rdname bf_priors
#' @export
bf_exponential <- function(rate = 1) {
  new_bf_prior("exponential", rate = rate)
}

#' @rdname bf_priors
#' @export
bf_lognormal <- function(meanlog = 0, sdlog = 1) {
  new_bf_prior("lognormal", meanlog = meanlog, sdlog = sdlog)
}

#' @rdname bf_priors
#' @export
bf_uniform <- function(lower = 0, upper = 1) {
  if (!(upper > lower)) cli::cli_abort("{.arg upper} must be greater than {.arg lower}.")
  new_bf_prior("uniform", lower = lower, upper = upper)
}

#' Empirical Heterogeneity Priors
#'
#' @description
#' Informative priors for \eqn{\tau} from large analyses of Cochrane reviews.
#' They are most useful when there are few studies, where \eqn{\tau} is
#' poorly identified by the data.
#'
#' \code{bf_turner()} gives the log-normal priors of Turner et al. (2015) for
#' log odds ratios (binary outcomes). \code{bf_rhodes()} gives the log-t
#' priors of Rhodes et al. (2015) for standardised mean differences
#' (continuous outcomes). The parameters come from
#' \code{bayesmeta::TurnerEtAlPrior()} and \code{bayesmeta::RhodesEtAlPrior()},
#' so the bayesmeta package must be installed. See those functions for the
#' outcome and comparison categories.
#'
#' @param outcome Outcome type, e.g. \code{"all-cause mortality"} for
#'   \code{bf_turner()} or \code{"mental health outcome"} for
#'   \code{bf_rhodes()}. \code{NA} gives the prior for a general outcome.
#' @param comparator1,comparator2 Intervention types compared: any of
#'   \code{"pharmacological"}, \code{"non-pharmacological"} and
#'   \code{"placebo / control"}.
#' @param area Medical area for \code{bf_rhodes()}: \code{"other"},
#'   \code{"respiratory"} or \code{"cancer"}.
#'
#' @return A \code{bf_prior} for \eqn{\tau}. See \code{\link{bf_priors}} for
#'   which backends support it.
#'
#' @references
#' Turner RM, Jackson D, Wei Y, Thompson SG, Higgins JPT (2015). Predictive
#' distributions for between-study heterogeneity and simple methods for their
#' application in Bayesian meta-analysis. \emph{Statistics in Medicine},
#' 34(6), 984--998.
#'
#' Rhodes KM, Turner RM, Higgins JPT (2015). Predictive distributions were
#' developed for the extent of heterogeneity in meta-analyses of continuous
#' outcome data. \emph{Journal of Clinical Epidemiology}, 68(1), 52--60.
#'
#' @examples
#' \dontrun{
#' bf_turner("all-cause mortality", "pharmacological", "placebo / control")
#' bf_rhodes("mental health outcome", "non-pharmacological", "placebo / control")
#' }
#'
#' @keywords internal
#' @export
bf_turner <- function(outcome = NA,
                      comparator1 = "pharmacological",
                      comparator2 = "placebo / control") {
  rlang::check_installed("bayesmeta", reason = "for the Turner et al. priors.")
  tp <- bayesmeta::TurnerEtAlPrior(outcome, comparator1, comparator2)
  # parameters are log-normal (mu, sigma) for tau^2 and tau; use the tau row
  prior <- bf_lognormal(meanlog = tp$parameters["tau", "mu"],
                        sdlog = tp$parameters["tau", "sigma"])
  prior$source <- "Turner"
  prior
}

#' @rdname bf_turner
#' @export
bf_rhodes <- function(outcome = NA,
                      comparator1 = "pharmacological",
                      comparator2 = "placebo / control",
                      area = "other") {
  rlang::check_installed("bayesmeta", reason = "for the Rhodes et al. priors.")
  rp <- bayesmeta::RhodesEtAlPrior(outcome, comparator1, comparator2, area)
  # parameters are log-t(5) (location, scale) for tau^2 and tau
  prior <- new_bf_prior("log_t", df = 5,
                        location = rp$parameters["tau", "location"],
                        scale = rp$parameters["tau", "scale"])
  prior$source <- "Rhodes"
  prior
}

#' @export
print.bf_prior <- function(x, ...) {
  cat(format_bf_prior(x), "\n")
}


# =========================================================================
# Families, parameters and backend support
# =========================================================================

#' @keywords internal
#' @noRd
new_bf_prior <- function(family, ...) {
  structure(list(family = family, ...), class = "bf_prior")
}

# Which parameter each family is a prior for
.bf_mu_families  <- c("normal", "student_t")
.bf_tau_families <- c("half_normal", "half_cauchy", "half_t", "exponential",
                      "lognormal", "log_t", "uniform")

# Families each backend can fit exactly
.bf_backend_support <- list(
  brms      = list(mu = .bf_mu_families, tau = .bf_tau_families),
  bayesmeta = list(mu = "normal", tau = .bf_tau_families),
  MetaStan  = list(mu = "normal", tau = c("half_normal", "half_cauchy", "uniform"))
)

#' Check a bf_prior suits the parameter and the backend
#'
#' @param prior A \code{bf_prior}, or NULL.
#' @param param "mu" or "tau".
#' @param backend Name in \code{.bf_backend_support}, or NULL to check the
#'   parameter only.
#'
#' @keywords internal
#' @noRd
check_bf_prior <- function(prior, param = c("mu", "tau"), backend = NULL) {
  param <- match.arg(param)
  if (is.null(prior)) return(invisible(prior))
  if (!inherits(prior, "bf_prior")) {
    cli::cli_abort(
      "{.arg {param}_prior} must be a {.cls bf_prior}, e.g. from {.fn bf_normal}."
    )
  }
  fam <- prior$family
  valid <- if (param == "mu") .bf_mu_families else .bf_tau_families
  if (!fam %in% valid) {
    other <- if (param == "mu") "tau" else "mu"
    suggest <- if (param == "mu") "bf_normal()" else "bf_half_normal()"
    cli::cli_abort(c(
      "A {.val {fam}} prior is for {other}, not {param}.",
      "i" = "For {param}, use e.g. {.code {suggest}}."
    ))
  }
  if (!is.null(backend)) {
    ok <- .bf_backend_support[[backend]][[param]]
    if (!fam %in% ok) {
      cli::cli_abort(c(
        "{backend} cannot fit a {.val {fam}} prior for {param}.",
        "i" = "Priors {backend} supports for {param}: {.val {ok}}.",
        "i" = "See {.help bayesfoRest::bf_priors} for the full table."
      ))
    }
    if (backend == "MetaStan" && param == "tau") {
      shift <- prior$mean %||% prior$location %||% prior$lower %||% 0
      if (shift != 0) {
        cli::cli_abort(
          "MetaStan tau priors must start at 0, so the mean, location or lower bound must be 0."
        )
      }
    }
  }
  if (fam == "uniform" && param == "tau" && prior$lower != 0) {
    cli::cli_abort("A uniform prior for tau must have {.code lower = 0}.")
  }
  invisible(prior)
}


# =========================================================================
# Densities (for plotting priors)
# =========================================================================

#' Density function for a bf_prior on the model scale
#'
#' Half distributions are truncated at 0 and renormalised.
#'
#' @keywords internal
#' @noRd
bf_prior_density <- function(prior) {
  switch(
    prior$family,
    normal      = function(x) stats::dnorm(x, prior$mean, prior$sd),
    student_t   = function(x) stats::dt((x - prior$location) / prior$scale,
                                        prior$df) / prior$scale,
    half_normal = truncate_at_zero(function(x) stats::dnorm(x, prior$mean, prior$sd)),
    half_cauchy = truncate_at_zero(function(x) stats::dcauchy(x, prior$location, prior$scale)),
    half_t      = truncate_at_zero(function(x) stats::dt((x - prior$location) / prior$scale,
                                                         prior$df) / prior$scale),
    exponential = function(x) stats::dexp(x, prior$rate),
    lognormal   = function(x) stats::dlnorm(x, prior$meanlog, prior$sdlog),
    log_t       = function(x) ifelse(
      x > 0,
      stats::dt((log(pmax(x, 1e-300)) - prior$location) / prior$scale,
                prior$df) / (prior$scale * pmax(x, 1e-300)),
      0),
    uniform     = function(x) stats::dunif(x, prior$lower, prior$upper),
    cli::cli_abort("Unsupported prior family {.val {prior$family}}.")
  )
}

#' Truncate a density at zero and renormalise it
#'
#' @keywords internal
#' @noRd
truncate_at_zero <- function(f) {
  mass <- stats::integrate(f, 0, Inf)$value
  function(x) ifelse(x < 0, 0, f(x) / mass)
}


# =========================================================================
# Translation to backend-specific formats
# =========================================================================

#' Translate bf_prior to a brms prior
#'
#' @param prior A \code{bf_prior} object.
#' @param class Character. brms prior class ("Intercept" for mu, "sd" for tau).
#'
#' @return A \code{brmsprior} object.
#'
#' @keywords internal
#' @noRd
bf_prior_to_brms <- function(prior, class = "Intercept") {
  rlang::check_installed("brms")

  # sd parameters are lower-bounded at 0 in brms, so half_* priors are the
  # untruncated family on class "sd"
  str <- switch(
    prior$family,
    normal      = ,
    half_normal = paste0("normal(", prior$mean, ", ", prior$sd, ")"),
    half_cauchy = paste0("cauchy(", prior$location, ", ", prior$scale, ")"),
    student_t   = ,
    half_t      = paste0("student_t(", prior$df, ", ", prior$location, ", ",
                         prior$scale, ")"),
    exponential = paste0("exponential(", prior$rate, ")"),
    lognormal   = paste0("lognormal(", prior$meanlog, ", ", prior$sdlog, ")"),
    uniform     = return(brms::prior_string(
      paste0("uniform(", prior$lower, ", ", prior$upper, ")"),
      class = class, ub = prior$upper)),
    # No log-t in Stan: flat sd prior plus the log-t density on log(sd) with
    # the Jacobian. Meta-analysis models have one sd, sd_1[1].
    log_t       = return(rbind(
      brms::prior_string("", class = "sd"),
      brms::prior_string(paste0(
        "target += student_t_lpdf(log(sd_1[1]) | ", prior$df, ", ",
        prior$location, ", ", prior$scale, ") - log(sd_1[1])"),
        check = FALSE))),
    cli::cli_abort("Unknown prior family {.val {prior$family}}.")
  )

  brms::prior_string(str, class = class)
}

#' Build Backend-Specific Prior Set from bf_prior Specifications
#'
#' @description
#' Translates a list containing \code{mu_prior} and \code{tau_prior}
#' (both \code{bf_prior} objects) into the format expected by a specific
#' model backend. Stops if the backend cannot fit a prior exactly.
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
  backend <- prior_backend(model)

  check_bf_prior(mu_prior, "mu", backend)
  check_bf_prior(tau_prior, "tau", backend)

  if (identical(backend, "brms")) {
    if (identical(tau_prior$family, "log_t") &&
        sum(grepl("^sd_", brms::variables(model))) != 1) {
      cli::cli_abort(
        "A log-t tau prior needs a brms model with a single study-level SD."
      )
    }
    rbind(
      if (!is.null(mu_prior)) bf_prior_to_brms(mu_prior, class = "Intercept"),
      if (!is.null(tau_prior)) bf_prior_to_brms(tau_prior, class = "sd")
    )

  } else if (identical(backend, "bayesmeta")) {
    list(
      mu_prior  = if (!is.null(mu_prior)) c(mu_prior$mean, mu_prior$sd),
      tau_prior = if (!is.null(tau_prior)) bf_prior_density(tau_prior)
    )

  } else if (identical(backend, "MetaStan")) {
    list(
      mu_prior       = if (!is.null(mu_prior)) c(mu_prior$mean, mu_prior$sd),
      # MetaStan's tau prior is a single scale (or upper bound) parameter
      tau_prior      = if (!is.null(tau_prior)) switch(
        tau_prior$family,
        half_normal = tau_prior$sd,
        half_cauchy = tau_prior$scale,
        uniform     = tau_prior$upper
      ),
      tau_prior_dist = if (!is.null(tau_prior)) switch(
        tau_prior$family,
        half_normal = "half-normal",
        half_cauchy = "half-cauchy",
        uniform     = "uniform"
      )
    )

  } else {
    cli::cli_abort(c(
      "Priors cannot be changed for {.cls {class(model)[1]}} models.",
      "i" = "Prior sensitivity refits support brms, bayesmeta and MetaStan models."
    ))
  }
}

#' Backend name used in .bf_backend_support, or NULL if priors can't be set
#'
#' @keywords internal
#' @noRd
prior_backend <- function(model) {
  if (inherits(model, "brmsfit")) "brms"
  else if (inherits(model, "bayesmeta")) "bayesmeta"
  else if (inherits(model, "meta_stan")) "MetaStan"
  else NULL
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

  # Unicode math-script letters (escaped to keep this file ASCII-only)
  script_N <- "\U0001D4A9"  # script N
  script_C <- "\U0001D49E"  # script C
  script_t <- "\U0001D4C9"  # script t
  script_E <- "\U2130"      # script E
  script_U <- "\U0001D4B0"  # script U
  superscript_plus <- "\U207A"  # superscript plus

  # Trailing thin space (U+2009) prevents STIX Two Math font from
  # rendering a doubled closing parenthesis in gt tables.
  thin <- "\U2009"

  # Subscript degrees of freedom for single-digit df
  sub_df <- function(df) {
    sub_digits <- c(
      "0" = "\U2080", "1" = "\U2081", "2" = "\U2082", "3" = "\U2083",
      "4" = "\U2084", "5" = "\U2085", "6" = "\U2086", "7" = "\U2087",
      "8" = "\U2088", "9" = "\U2089"
    )
    if (grepl("^[0-9]$", df)) sub_digits[[as.character(df)]] else as.character(df)
  }
  num <- function(x) format(signif(x, 3))

  label <- switch(
    prior$family,

    normal = paste0(script_N, "(", prior$mean, ",", prior$sd, ")"),

    student_t = paste0(script_t, sub_df(prior$df),
                       "(", prior$location, ",", prior$scale, ")"),

    half_normal = paste0(script_N, superscript_plus,
                         "(", prior$mean, ",", prior$sd, ")"),

    half_cauchy = paste0(script_C, superscript_plus,
                         "(", prior$location, ",", prior$scale, ")"),

    half_t = paste0(script_t, sub_df(prior$df), superscript_plus,
                    "(", prior$location, ",", prior$scale, ")"),

    exponential = paste0(script_E, "(", prior$rate, ")"),

    lognormal = paste0("log", script_N, "(", num(prior$meanlog), ",",
                       num(prior$sdlog), ")"),

    log_t = paste0("log", script_t, sub_df(prior$df), "(", num(prior$location),
                   ",", num(prior$scale), ")"),

    uniform = paste0(script_U, "(", prior$lower, ",", prior$upper, ")"),

    paste0("Unknown(", prior$family, ")")
  )

  if (!is.null(prior$source)) label <- paste0(prior$source, " ", label)
  paste0(label, thin)
}
