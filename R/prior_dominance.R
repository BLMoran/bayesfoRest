# =============================================================================
# prior_dominance(): quantify how much a Bayesian meta-analysis is driven by the
# prior vs the data -- backend-agnostic.
#
# Priors are obtained as DRAWS (from the fit itself, or supplied directly), so
# no distribution needs to be re-specified and any modelling package is
# supported. Everything downstream is distribution-free:
#
#   contraction = 1 - Var(posterior)/Var(prior)   (~1 data-dominated, ~0 prior-dominated)
#   percentile  = empirical prior CDF at the posterior median
#                 (interior = no conflict; near 0/1 = possible prior-data conflict)
#   tau grid    = pooled mu (+95% CrI) across a set of tau-prior scales
#                 (brms / rstanarm refits; skipped for other backends)
#
# Prior draws:
#   * pass `prior_draws` (a data frame / draws object with the mu & tau columns)
#     -- works for ANY backend (brms, cmdstanr, rstan, JAGS, ...); or
#   * let the function sample them from the fit (brms: sample_prior = "only";
#     rstanarm: prior_PD = TRUE).
# =============================================================================

#' Assess prior dominance in a Bayesian meta-analysis (backend-agnostic)
#'
#' @param fit A fitted model. Posterior draws are read via
#'   \code{posterior::as_draws_df()}, so any backend that supports it works.
#' @param prior_draws Optional PRIOR draws (data frame / draws object) with the
#'   mu and tau columns. If supplied it is used directly (portable path for
#'   cmdstanr / rstan / JAGS). If NULL, priors are sampled from the fit.
#' @param mu_par,tau_par Column names for the pooled effect and between-study SD.
#'   Auto-detected when NULL.
#' @param tau_grid Half-normal scales for the tau-prior sensitivity refits
#'   (brms / rstanarm only). NULL to skip.
#' @param refresh Passed to refits.
#'
#' @return An object of class \code{prior_dominance}.
#' @keywords internal
#' @export
#'
#' @importFrom stats var median quantile update
prior_dominance <- function(fit,
                            prior_draws = NULL,
                            mu_par = NULL,
                            tau_par = NULL,
                            tau_grid = c(0.5, 0.75, 1.5),
                            refresh = 0) {

  post <- posterior::as_draws_df(fit)

  mu_par  <- mu_par  %||% intersect(c("b_Intercept", "mu", "pooled"), names(post))[1]
  tau_par <- tau_par %||% c(grep("^sd_.*__Intercept$", names(post), value = TRUE),
                            intersect("tau", names(post)))[1]
  if (is.na(mu_par) || is.na(tau_par))
    stop("Could not locate mu / tau draws. Pass mu_par / tau_par explicitly.")

  post_mu  <- post[[mu_par]]
  post_tau <- post[[tau_par]]

  # ---- prior draws: supplied, or sampled from the fit --------------------
  if (is.null(prior_draws))
    prior_draws <- .pd_prior_draws(fit, c(mu_par, tau_par), refresh = refresh)
  prior_draws <- posterior::as_draws_df(prior_draws)
  if (!all(c(mu_par, tau_par) %in% names(prior_draws)))
    stop("`prior_draws` must contain columns '", mu_par, "' and '", tau_par, "'.")
  prior_mu  <- prior_draws[[mu_par]]
  prior_tau <- prior_draws[[tau_par]]

  # ---- 1. contraction (guarding improper / infinite-variance priors) -----
  contraction_of <- function(prior, post) {
    pv <- stats::var(prior, na.rm = TRUE)
    if (!is.finite(pv) || pv == 0) return(NA_real_)
    1 - stats::var(post, na.rm = TRUE) / pv
  }
  interp <- function(c) dplyr::case_when(
    is.na(c)  ~ "improper prior - undefined",
    c >= 0.75 ~ "data-dominated",
    c >= 0.30 ~ "informed by both",
    TRUE      ~ "prior-dominated")

  cmu  <- contraction_of(prior_mu,  post_mu)
  ctau <- contraction_of(prior_tau, post_tau)
  contraction <- tibble::tibble(
    parameter      = c("mu", "tau"),
    prior_var      = c(stats::var(prior_mu), stats::var(prior_tau)),
    posterior_var  = c(stats::var(post_mu),  stats::var(post_tau)),
    contraction    = c(cmu, ctau),
    interpretation = interp(c(cmu, ctau)))

  # ---- 2. empirical prior percentile at the posterior median -------------
  pct  <- function(prior, x) mean(prior <= x, na.rm = TRUE)
  pmu  <- pct(prior_mu,  stats::median(post_mu))
  ptau <- pct(prior_tau, stats::median(post_tau))
  percentile <- tibble::tibble(
    parameter        = c("mu", "tau"),
    posterior_median = c(stats::median(post_mu), stats::median(post_tau)),
    prior_percentile = c(pmu, ptau),
    flag             = dplyr::case_when(
      c(pmu, ptau) > 0.95 ~ "upper tail - possible prior-data conflict",
      c(pmu, ptau) < 0.05 ~ "lower tail - possible prior-data conflict",
      TRUE                ~ "interior - no conflict"))

  # ---- 3. tau-prior sensitivity (brms / rstanarm refits) -----------------
  tau_sensitivity <- NULL
  if (!is.null(tau_grid) && (inherits(fit, "brmsfit") || inherits(fit, "stanreg")))
    tau_sensitivity <- .pd_tau_grid(fit, mu_par, tau_grid, refresh)

  # ---- verdict -----------------------------------------------------------
  mu_range <- if (!is.null(tau_sensitivity)) diff(range(tau_sensitivity$mu_median)) else NA_real_
  verdict <- sprintf(
    "mu contraction = %s (%s); tau contraction = %s (%s).%s",
    ifelse(is.na(cmu), "NA", sprintf("%.2f", cmu)),  interp(cmu),
    ifelse(is.na(ctau), "NA", sprintf("%.2f", ctau)), interp(ctau),
    if (!is.na(mu_range)) sprintf(" Pooled mu shifts by %.3f across the tau-prior grid.", mu_range) else "")

  # flag any prior-data conflict detected in the tau grid
  if (!is.null(tau_sensitivity) && any(tau_sensitivity$conflict != "compatible")) {
    bad <- tau_sensitivity$tau_prior[tau_sensitivity$conflict != "compatible"]
    verdict <- paste0(verdict,
                      sprintf(" Prior-data conflict flagged for: %s (data want more heterogeneity than the prior allows).",
                              paste(bad, collapse = ", ")))
  } else if (!is.null(tau_sensitivity)) {
    verdict <- paste0(verdict, " No prior-data conflict across the tau-prior grid.")
  }

  structure(list(contraction = contraction, percentile = percentile,
                 tau_sensitivity = tau_sensitivity, verdict = verdict,
                 mu_par = mu_par, tau_par = tau_par),
            class = "prior_dominance")
}


# ---------------------------------------------------------------------------
# Prior draws from the fit (dispatch); prior_draws= bypasses all of this.
# ---------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.pd_prior_draws <- function(fit, pars, refresh = 0) UseMethod(".pd_prior_draws")

#' @keywords internal
#' @noRd
.pd_prior_draws.brmsfit <- function(fit, pars, refresh = 0) {
  fp <- stats::update(fit, sample_prior = "only", refresh = refresh, silent = 2)
  posterior::as_draws_df(fp)
}

#' @keywords internal
#' @noRd
.pd_prior_draws.stanreg <- function(fit, pars, refresh = 0) {
  fp <- stats::update(fit, prior_PD = TRUE, refresh = refresh)
  posterior::as_draws_df(as.data.frame(fp))
}

#' @keywords internal
#' @noRd
.pd_prior_draws.default <- function(fit, pars, refresh = 0)
  stop("Automatic prior sampling isn't available for this model class.\n",
       "Supply prior draws via `prior_draws = <data.frame with ",
       paste(pars, collapse = ', '), ">`.")


# ---------------------------------------------------------------------------
# tau-prior grid (brms / rstanarm)
# ---------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.pd_tau_grid <- function(fit, mu_par, tau_grid, refresh) UseMethod(".pd_tau_grid")

#' @keywords internal
#' @noRd
.pd_tau_grid.brmsfit <- function(fit, mu_par, tau_grid, refresh) {
  d0  <- posterior::as_draws_df(fit)
  grp <- grep("^sd_.*__Intercept$", names(d0), value = TRUE)
  group_name <- if (length(grp)) sub("^sd_(.*)__Intercept$", "\\1", grp[1]) else "Author"
  mu_prior <- brms::prior_summary(fit)
  mu_row   <- mu_prior[mu_prior$class == "Intercept", "prior"]
  mu_str   <- if (length(mu_row) && nzchar(mu_row[1])) mu_row[1] else "normal(0, 1)"

  purrr::map_dfr(tau_grid, function(s) {
    pr <- c(brms::prior_string(mu_str, class = "Intercept"),
            brms::prior_string(sprintf("normal(0, %g)", s),
                               class = "sd", group = group_name, lb = "0"))
    f  <- stats::update(fit, prior = pr, refresh = refresh, silent = 2)
    dd <- posterior::as_draws_df(f)
    mm <- dd[[mu_par]]
    tt <- dd[[grep("^sd_.*__Intercept$", names(dd), value = TRUE)[1]]]

    tau_med <- stats::median(tt)
    # where does posterior tau sit in THIS row's HN(0, s) prior?
    # folded-normal CDF at tau_med (tau >= 0). Upper tail => prior too tight.
    tau_pctile <- 2 * stats::pnorm(tau_med / s) - 1
    # For a variance parameter, only the UPPER tail is a genuine conflict:
    # the prior is squeezing tau toward zero against the data's wishes.
    conflict <- dplyr::if_else(
      tau_pctile > 0.95,
      "conflict: prior too tight (data want more heterogeneity)",
      "compatible")

    tibble::tibble(tau_prior = sprintf("HN(0, %g)", s),
                   mu_median = stats::median(mm),
                   mu_lo = stats::quantile(mm, .025, names = FALSE),
                   mu_hi = stats::quantile(mm, .975, names = FALSE),
                   pr_benefit = mean(mm < 0),
                   tau_median = tau_med,
                   tau_prior_pctile = tau_pctile,
                   conflict = conflict)
  })
}

#' @keywords internal
#' @noRd
.pd_tau_grid.stanreg <- function(fit, mu_par, tau_grid, refresh) {
  warning("tau-prior grid not implemented for stanreg; skipping.")
  NULL
}


# ---------------------------------------------------------------------------
# print + gt
# ---------------------------------------------------------------------------
#' @export
print.prior_dominance <- function(x, ...) {
  cat("Prior dominance assessment\n==========================\n\n")
  cat("Posterior contraction (1 = data-dominated, 0 = prior-dominated):\n")
  print(x$contraction); cat("\n")
  cat("Posterior median within its prior:\n")
  print(x$percentile); cat("\n")
  if (!is.null(x$tau_sensitivity)) {
    cat("Pooled effect across tau-prior grid:\n")
    print(x$tau_sensitivity); cat("\n")
  }
  cat("Verdict:\n ", x$verdict, "\n")
  invisible(x)
}

#' Render a prior_dominance object as a single gt table
#' @param x A \code{prior_dominance} object.
#' @param decimals Decimal places for numeric columns.
#' @return A \code{gt_tbl}.
#' @keywords internal
#' @export
gt_prior_dominance <- function(x, decimals = 3) {

  contr <- x$contraction |>
    dplyr::transmute(
      Section = "Contraction",
      Item    = parameter,
      Value   = ifelse(is.na(contraction), "\u2014",
                       sprintf(paste0("%.", decimals, "f"), contraction)),
      Detail  = interpretation)

  perc <- x$percentile |>
    dplyr::transmute(
      Section = "Prior percentile of posterior median",
      Item    = parameter,
      Value   = sprintf("%.1f%%", 100 * prior_percentile),
      Detail  = flag)

  tau <- if (!is.null(x$tau_sensitivity)) {
    x$tau_sensitivity |>
      dplyr::transmute(
        Section = "Pooled effect across tau priors",
        Item    = tau_prior,
        Value   = sprintf(paste0("%.", decimals, "f [%.", decimals, "f, %.", decimals, "f]"),
                          mu_median, mu_lo, mu_hi),
        Detail  = sprintf("P(benefit) = %.2f; tau = %.2f (%.0f%% of prior; %s)",
                          pr_benefit, tau_median, 100 * tau_prior_pctile, conflict))
  } else NULL

  dplyr::bind_rows(contr, perc, tau) |>
    gt::gt(groupname_col = "Section", rowname_col = "Item") |>
    gt::tab_header(title = gt::md("**Prior dominance assessment**"),
                   subtitle = x$verdict) |>
    gt::cols_label(Value = "Estimate", Detail = "") |>
    gt::tab_options(row_group.font.weight = "bold") |>
    gt::sub_missing(missing_text = "\u2014")
}

`%||%` <- function(a, b) if (is.null(a)) b else a


# =============================================================================
# Continuous prior-scale sweep (the "power-prior" idea applied to the mu prior):
# tighten/loosen the mu prior SD and watch whether the posterior tracks the data
# or gets dragged toward the prior mean. A flat mu-vs-prior_sd curve = the data
# hold the estimate; a curve bending toward the prior mean as the prior tightens
# = prior dominance made visible.
# =============================================================================

#' Continuous mu-prior scale sweep
#'
#' @param fit A brms / rstanarm meta-analysis fit.
#' @param scales Numeric vector of mu-prior SDs to sweep.
#' @param prior_mean Mean of the mu prior held fixed across the sweep
#'   (default 0). Set to your informative-prior mean to test how hard a tight
#'   informative prior pulls the estimate.
#' @param mu_par Pooled-effect draw name (auto-detected).
#' @param refresh Passed to refits.
#' @return A tibble: prior_sd, mu_median, mu_lo, mu_hi, pr_benefit.
#' @keywords internal
#' @export
prior_scale_sweep <- function(fit,
                              scales = c(0.1, 0.25, 0.5, 1, 2, 5),
                              prior_mean = 0,
                              mu_par = NULL,
                              refresh = 0) {
  stopifnot(inherits(fit, "brmsfit") || inherits(fit, "stanreg"))
  d0 <- posterior::as_draws_df(fit)
  mu_par <- mu_par %||% intersect(c("b_Intercept", "mu", "pooled"), names(d0))[1]

  # keep the model's tau prior unchanged; only vary the mu prior
  tau_rows <- brms::prior_summary(fit)
  tau_str  <- tau_rows[tau_rows$class == "sd", "prior"]
  tau_str  <- if (length(tau_str) && nzchar(tau_str[1])) tau_str[1] else "normal(0, 0.75)"
  grp      <- grep("^sd_.*__Intercept$", names(d0), value = TRUE)
  group_nm <- if (length(grp)) sub("^sd_(.*)__Intercept$", "\\1", grp[1]) else "Author"

  purrr::map_dfr(scales, function(s) {
    pr <- c(brms::prior_string(sprintf("normal(%g, %g)", prior_mean, s), class = "Intercept"),
            brms::prior_string(tau_str, class = "sd", group = group_nm, lb = "0"))
    f  <- stats::update(fit, prior = pr, refresh = refresh, silent = 2)
    mm <- posterior::as_draws_df(f)[[mu_par]]
    tibble::tibble(prior_sd = s,
                   mu_median = stats::median(mm),
                   mu_lo = stats::quantile(mm, .025, names = FALSE),
                   mu_hi = stats::quantile(mm, .975, names = FALSE),
                   pr_benefit = mean(mm < 0))
  })
}

#' Plot a prior-scale sweep (mu vs prior SD)
#' @param sweep Output of \code{prior_scale_sweep()}.
#' @param prior_mean The prior mean used (drawn as a reference line).
#' @param x_label,y_label Axis labels.
#' @param font Font family.
#' @return A ggplot.
#' @keywords internal
#' @export
plot_prior_scale_sweep <- function(sweep, prior_mean = 0,
                                   x_label = "mu-prior SD (log scale)",
                                   y_label = "Pooled effect",
                                   font = NULL) {
  ggplot2::ggplot(sweep, ggplot2::aes(x = prior_sd, y = mu_median)) +
    ggplot2::geom_hline(yintercept = prior_mean, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = mu_lo, ymax = mu_hi),
                         fill = "dodgerblue", alpha = 0.15) +
    ggplot2::geom_line(colour = "blue", linewidth = 0.8) +
    ggplot2::geom_point(colour = "blue", size = 2) +
    ggplot2::scale_x_log10() +
    ggplot2::annotate("text", x = min(sweep$prior_sd), y = prior_mean,
                      label = "prior mean", hjust = 0, vjust = -0.5,
                      size = 3, colour = "grey40", family = font) +
    ggplot2::labs(x = x_label, y = y_label,
                  subtitle = "Flat curve => estimate is data-driven; bending toward the prior mean as the prior tightens => prior influence") +
    ggplot2::theme_light() +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 9, family = font),
                   axis.title = ggplot2::element_text(family = font),
                   axis.text = ggplot2::element_text(family = font))
}


# =============================================================================
# Prior effective sample size (ESS), in study-equivalents.
# For a mean-like parameter, a prior with variance v_prior carries precision
# 1/v_prior; expressed relative to the typical within-study precision it is
# "worth" that many studies. This makes prior informativeness concrete:
# "the prior contributes the equivalent of N studies" against your real k.
# =============================================================================

#' Prior effective sample size for the pooled effect (study-equivalents)
#'
#' @param fit A fitted meta-analysis model.
#' @param data The study-level data (needs the sei column) to define a typical
#'   study precision. Alternatively pass \code{typical_se} directly.
#' @param typical_se Optional numeric: the representative within-study SE. If
#'   NULL, computed as the pooled (harmonic-mean-precision) SE from \code{data$sei}.
#' @param prior_draws,mu_par As in \code{prior_dominance()}.
#' @return A list: prior_sd, typical_se, prior_ess (study-equivalents), k (actual
#'   studies), and the prior's share of total precision.
#' @keywords internal
#' @export
prior_ess <- function(fit, data = NULL, typical_se = NULL,
                      prior_draws = NULL, mu_par = NULL) {

  post <- posterior::as_draws_df(fit)
  mu_par <- mu_par %||% intersect(c("b_Intercept", "mu", "pooled"), names(post))[1]

  if (is.null(prior_draws)) prior_draws <- .pd_prior_draws(fit, mu_par)
  prior_draws <- posterior::as_draws_df(prior_draws)
  prior_sd <- stats::sd(prior_draws[[mu_par]], na.rm = TRUE)

  if (is.null(typical_se)) {
    if (is.null(data) || is.null(data$sei))
      stop("Provide `data` (with a sei column) or `typical_se`.")
    # precision-equivalent single-study SE: 1/mean(1/sei^2) -> sqrt
    typical_se <- sqrt(1 / mean(1 / data$sei^2, na.rm = TRUE))
    k <- sum(is.finite(data$sei))
  } else k <- NA_integer_

  prior_prec  <- 1 / prior_sd^2
  study_prec  <- 1 / typical_se^2
  prior_ess   <- prior_prec / study_prec                  # study-equivalents
  data_prec   <- if (!is.na(k)) k * study_prec else NA_real_
  prior_share <- if (!is.na(data_prec)) prior_prec / (prior_prec + data_prec) else NA_real_

  structure(list(prior_sd = prior_sd, typical_se = typical_se,
                 prior_ess = prior_ess, k = k, prior_share = prior_share),
            class = "prior_ess")
}

#' @export
print.prior_ess <- function(x, ...) {
  cat(sprintf(
    "Prior effective sample size (pooled effect)\n  prior SD           = %.3f\n  typical study SE   = %.3f\n  prior ESS          = %.2f study-equivalents\n",
    x$prior_sd, x$typical_se, x$prior_ess))
  if (!is.na(x$k))
    cat(sprintf("  actual studies (k) = %d\n  prior share of total precision = %.1f%%\n",
                x$k, 100 * x$prior_share))
  invisible(x)
}


# =============================================================================
# Continuous TAU-prior scale sweep -- the tau analogue of prior_scale_sweep().
# Vary the half-normal scale on tau (holding the mu prior fixed) and track BOTH
# the pooled effect mu and the heterogeneity tau. Shows whether the pooled
# estimate and the heterogeneity estimate depend on how the tau prior is set.
# =============================================================================

#' Continuous tau-prior scale sweep
#'
#' @param fit A brms / rstanarm meta-analysis fit.
#' @param scales Numeric vector of half-normal scales for the tau prior.
#' @param mu_par Pooled-effect draw name (auto-detected).
#' @param refresh Passed to refits.
#' @return A tibble: tau_prior_sd, mu_median, mu_lo, mu_hi, pr_benefit,
#'   tau_median, tau_lo, tau_hi, tau_prior_pctile, conflict.
#' @keywords internal
#' @export
tau_scale_sweep <- function(fit,
                            scales = c(0.25, 0.5, 0.75, 1, 1.5, 3),
                            mu_par = NULL,
                            refresh = 0) {
  stopifnot(inherits(fit, "brmsfit") || inherits(fit, "stanreg"))
  d0 <- posterior::as_draws_df(fit)
  mu_par <- mu_par %||% intersect(c("b_Intercept", "mu", "pooled"), names(d0))[1]

  # hold the mu prior fixed at the fitted value
  mu_rows <- brms::prior_summary(fit)
  mu_str  <- mu_rows[mu_rows$class == "Intercept", "prior"]
  mu_str  <- if (length(mu_str) && nzchar(mu_str[1])) mu_str[1] else "normal(0, 1)"
  grp     <- grep("^sd_.*__Intercept$", names(d0), value = TRUE)
  group_nm <- if (length(grp)) sub("^sd_(.*)__Intercept$", "\\1", grp[1]) else "Author"

  purrr::map_dfr(scales, function(s) {
    pr <- c(brms::prior_string(mu_str, class = "Intercept"),
            brms::prior_string(sprintf("normal(0, %g)", s),
                               class = "sd", group = group_nm, lb = "0"))
    f  <- stats::update(fit, prior = pr, refresh = refresh, silent = 2)
    dd <- posterior::as_draws_df(f)
    mm <- dd[[mu_par]]
    tt <- dd[[grep("^sd_.*__Intercept$", names(dd), value = TRUE)[1]]]

    tau_med    <- stats::median(tt)
    tau_pctile <- 2 * stats::pnorm(tau_med / s) - 1          # folded-normal CDF
    conflict   <- dplyr::if_else(tau_pctile > 0.95,
                                 "conflict: prior too tight", "compatible")

    tibble::tibble(
      tau_prior_sd = s,
      mu_median = stats::median(mm),
      mu_lo = stats::quantile(mm, .025, names = FALSE),
      mu_hi = stats::quantile(mm, .975, names = FALSE),
      pr_benefit = mean(mm < 0),
      tau_median = tau_med,
      tau_lo = stats::quantile(tt, .025, names = FALSE),
      tau_hi = stats::quantile(tt, .975, names = FALSE),
      tau_prior_pctile = tau_pctile,
      conflict = conflict)
  })
}

#' Plot a tau-prior scale sweep (mu and tau vs tau-prior SD)
#'
#' @param sweep Output of \code{tau_scale_sweep()}.
#' @param mu_label,tau_label Facet labels for the two panels.
#' @param null_value Reference line for mu (default 0).
#' @param font Font family.
#' @return A ggplot (two stacked panels: mu and tau vs tau-prior SD).
#' @keywords internal
#' @export
plot_tau_scale_sweep <- function(sweep,
                                 mu_label = "Pooled effect (\u03bc)",
                                 tau_label = "Heterogeneity (\u03c4)",
                                 null_value = 0,
                                 font = NULL) {

  mu_df <- sweep |>
    dplyr::transmute(tau_prior_sd, panel = mu_label,
                     est = mu_median, lo = mu_lo, hi = mu_hi)
  tau_df <- sweep |>
    dplyr::transmute(tau_prior_sd, panel = tau_label,
                     est = tau_median, lo = tau_lo, hi = tau_hi)
  df <- dplyr::bind_rows(mu_df, tau_df) |>
    dplyr::mutate(panel = factor(panel, levels = c(mu_label, tau_label)))

  # null line only in the mu panel
  hlines <- data.frame(panel = factor(mu_label, levels = c(mu_label, tau_label)),
                       yint = null_value)

  ggplot2::ggplot(df, ggplot2::aes(x = tau_prior_sd, y = est)) +
    ggplot2::geom_hline(data = hlines, ggplot2::aes(yintercept = yint),
                        linetype = "dashed", colour = "grey50") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                         fill = "dodgerblue", alpha = 0.15) +
    ggplot2::geom_line(colour = "blue", linewidth = 0.8) +
    ggplot2::geom_point(colour = "blue", size = 2) +
    ggplot2::facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      x = "\u03c4-prior scale (half-normal SD, log scale)",
      y = NULL,
      subtitle = "Flat \u03bc across \u03c4-prior scales => pooled effect robust to the heterogeneity prior") +
    ggplot2::theme_light() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", family = font),
      plot.subtitle = ggplot2::element_text(size = 9, family = font),
      axis.title = ggplot2::element_text(family = font),
      axis.text = ggplot2::element_text(family = font),
      panel.grid.minor = ggplot2::element_blank())
}
