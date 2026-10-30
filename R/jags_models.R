#' Generate JAGS Model Files for Bayesian Meta-Analysis
#'
#' @description
#' Writes a JAGS model file for common meta-analysis models (binary, continuous,
#' count outcomes) with a user-specified prior on the between-study standard
#' deviation (sigma). These model files can be passed directly to
#' [R2jags::jags()].
#'
#' @param outcome Character: `"bin"`, `"ctns"`, or `"count"`.
#' @param model Character: model specification string.
#'   - Binary: `"std.norm"`, `"std.dt"`, `"reg.norm"`, `"reg.dt"`
#'   - Continuous: `"std.ta"`, `"std.mv"`, `"reg.ta"`, `"reg.mv"`
#'   - Count: `"std"`, `"std.unif"`, `"std.hc"`, `"reg"`, `"reg.unif"`,
#'     `"reg.hc"`
#' @param type Character: `"fix"` (common-effects) or `"ran"` (random-effects).
#' @param sigma_prior Character: a block of JAGS code defining `sigma` and
#'   `tau`. Must define both variables. Defaults to a half-normal(0, 0.5) prior.
#'   Use the helper functions:
#'   - `jags_prior_halfnormal(scale)` (default)
#'   - `jags_prior_halfcauchy(scale)`
#'   - `jags_prior_uniform(upper)`
#' @param path File path to write the model to. Defaults to a tempfile.
#'
#' @return The file path (invisibly), for passing to [R2jags::jags()].
#'
#' @keywords internal
#' @export
jags_model_file <- function(outcome = c("bin", "ctns", "count"),
                            model,
                            type = c("ran", "fix"),
                            sigma_prior = jags_prior_halfnormal(0.5),
                            path = NULL) {
  outcome <- rlang::arg_match(outcome)
  type <- rlang::arg_match(type)
  path <- path %||% tempfile(fileext = ".txt")

  likelihood <- .jags_likelihood(outcome, model, type)
  prior_mu <- .jags_mu_prior(outcome, model)
  derived <- .jags_derived(outcome)

  model_code <- paste0(
    "model {\n",
    likelihood, "\n",
    "  ## Priors\n",
    prior_mu, "\n",
    sigma_prior, "\n",
    derived, "\n",
    "}\n"
  )

  writeLines(model_code, path)
  invisible(path)
}


#' @describeIn jags_model_file Half-normal prior on sigma: `sigma ~ N(0, scale^2) T(0,)`
#' @param scale Scale parameter (SD of the half-normal). Default 0.5.
#' @export
jags_prior_halfnormal <- function(scale = 0.5) {
  tau_val <- 1 / scale^2
  paste0(
    "  sigma ~ dnorm(0, ", tau_val, ") T(0,)\n",
    "  tau <- pow(sigma, -2)\n"
  )
}

#' @describeIn jags_model_file Half-Cauchy prior on sigma via the
#'   scale-mixture representation.
#' @param scale Scale parameter of the half-Cauchy. Default 0.5.
#' @export
jags_prior_halfcauchy <- function(scale = 0.5) {
  paste0(
    "  sigma <- abs(z.sigma) / pow(epsilon.sigma, 0.5)\n",
    "  z.sigma ~ dnorm(0, tau.z.sigma)\n",
    "  epsilon.sigma ~ dgamma(0.5, 0.5)\n",
    "  tau.z.sigma <- pow(B.sigma, -2)\n",
    "  B.sigma <- ", scale, "\n",
    "  tau <- pow(sigma, -2)\n"
  )
}

#' @describeIn jags_model_file Uniform prior on sigma.
#' @param upper Upper bound of the uniform. Default 5.
#' @export
jags_prior_uniform <- function(upper = 5) {
  paste0(
    "  sigma ~ dunif(0, ", upper, ")\n",
    "  tau <- pow(sigma, -2)\n"
  )
}


# ── Internal: likelihood blocks ──────────────────────────────────────────────

.jags_likelihood <- function(outcome, model, type) {
  is_re <- identical(type, "ran")
  is_reg <- grepl("^reg", model)

  switch(outcome,
    bin = .jags_lik_bin(is_re, is_reg),
    ctns = .jags_lik_ctns(model, is_re, is_reg),
    count = .jags_lik_count(is_re, is_reg)
  )
}


.jags_lik_bin <- function(is_re, is_reg) {
  reg_block <- if (is_reg) {
    paste0(
      "    logit(pi0[s]) <- alpha[s] + Z0[s,] %*% beta0\n",
      "    logit(pi1[s]) <- alpha[s] + Z0[s,] %*% beta0 + delta[s]\n"
    )
  } else {
    paste0(
      "    logit(pi0[s]) <- alpha[s]\n",
      "    logit(pi1[s]) <- alpha[s] + delta[s]\n"
    )
  }

  delta_line <- if (is_re) {
    "    delta[s] ~ dnorm(mu, tau)\n"
  } else {
    "    delta[s] ~ dnorm(mu, 0.0001)\n"
  }

  paste0(
    "  for (s in 1:S) {\n",
    "    y0[s] ~ dbin(pi0[s], n0[s])\n",
    "    y1[s] ~ dbin(pi1[s], n1[s])\n",
    reg_block,
    "    alpha[s] ~ dnorm(0, 0.0001)\n",
    delta_line,
    "    gamma[s] <- exp(delta[s])\n",
    "  }\n"
  )
}


.jags_lik_ctns <- function(model, is_re, is_reg) {
  is_two_arm <- grepl("ta$", model)

  delta_line <- if (is_re) {
    if (!is_two_arm && is_reg) {
      "    alpha[s] ~ dnorm(mu, tau)\n"
    } else {
      "    delta[s] ~ dnorm(mu, tau)\n"
    }
  } else {
    "    delta[s] ~ dnorm(mu, 0.0001)\n"
  }

  if (is_two_arm) {
    reg_block <- if (is_reg) {
      paste0(
        "    alpha0[s] <- gamma0[s] + Z0[s,] %*% beta0\n",
        "    gamma0[s] ~ dnorm(0, 0.0001)\n"
      )
    } else {
      "    alpha0[s] ~ dnorm(0, 0.0001)\n"
    }

    paste0(
      "  for (s in 1:S) {\n",
      "    y0[s] ~ dnorm(alpha0[s], prec0[s])\n",
      "    y1[s] ~ dnorm(alpha1[s], prec1[s])\n",
      "    alpha1[s] <- alpha0[s] + delta[s]\n",
      "    prec0[s] <- pow(se0[s], -2)\n",
      "    prec1[s] <- pow(se1[s], -2)\n",
      reg_block,
      delta_line,
      "  }\n"
    )
  } else {
    if (is_reg) {
      paste0(
        "  for (s in 1:S) {\n",
        "    y[s] ~ dnorm(delta[s], prec[s])\n",
        "    delta[s] <- alpha[s] + Z0[s,] %*% beta0\n",
        delta_line,
        "  }\n"
      )
    } else {
      paste0(
        "  for (s in 1:S) {\n",
        "    y[s] ~ dnorm(delta[s], prec[s])\n",
        delta_line,
        "  }\n"
      )
    }
  }
}


.jags_lik_count <- function(is_re, is_reg) {
  reg_block <- if (is_reg) {
    paste0(
      "    log(lambda0[s]) <- xi0[s] + Z0[s,] %*% beta0 + log(p0[s])\n",
      "    log(lambda1[s]) <- xi1[s] + log(p1[s])\n"
    )
  } else {
    paste0(
      "    log(lambda0[s]) <- xi0[s] + log(p0[s])\n",
      "    log(lambda1[s]) <- xi1[s] + log(p1[s])\n"
    )
  }

  delta_line <- if (is_re) {
    "    delta[s] ~ dnorm(mu, tau)\n"
  } else {
    "    delta[s] ~ dnorm(mu, 0.0001)\n"
  }

  paste0(
    "  for (s in 1:S) {\n",
    "    y0[s] ~ dpois(lambda0[s])\n",
    "    y1[s] ~ dpois(lambda1[s])\n",
    reg_block,
    "    xi0[s] ~ dunif(-5, 5)\n",
    "    xi1[s] <- xi0[s] + delta[s]\n",
    delta_line,
    "    gamma[s] <- exp(delta[s])\n",
    "  }\n"
  )
}


# ── Internal: mu prior and derived quantities ────────────────────────────────

.jags_mu_prior <- function(outcome, model) {
  is_reg <- grepl("^reg", model)
  is_dt <- grepl("dt$", model)

  reg_line <- if (is_reg) {
    "  beta0[1:J] ~ dnorm(m.beta0[], tau.beta0[,])\n"
  } else {
    ""
  }

  mu_line <- if (is_dt) {
    paste0(
      "  mu ~ dt(0, 0.5, v)\n",
      "  v ~ dunif(0, 8)\n"
    )
  } else {
    "  mu ~ dnorm(0, 0.0001)\n"
  }

  paste0(reg_line, mu_line)
}


.jags_derived <- function(outcome) {
  switch(outcome,
    bin   = "  rho <- exp(mu)\n",
    count = "  IRR <- exp(mu)\n",
    ctns  = ""
  )
}
