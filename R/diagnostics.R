# =============================================================================
# Model-agnostic MCMC diagnostics for bayesfoRest.
#
# One generic, dispatching on model class the same way bayes_forest() does.
# Works on any fit bayesfoRest supports: brmsfit, stanfit (rstan / MetaStan),
# CmdStanMCMC (cmdstanr), stanreg. Builds a single-page dashboard directly from
# bayesplot / posterior / ggplot2 -- no external panel helpers, no glyph fonts.
#
# Display features:
#   * pooled effect / heterogeneity relabelled as mu / tau (per backend),
#     with a user override via `labels`.
#   * degenerate residual sigma (fixed by se()) excluded from trace/ACF.
#   * pass/fail table encoded by cell fill + plain words (device-independent).
# =============================================================================

#' Single-Page MCMC Diagnostics
#'
#' @param object A fitted model: brmsfit, stanfit (rstan / MetaStan),
#'   CmdStanMCMC (cmdstanr), or stanreg.
#' @param pars Character vector of parameters for the trace / autocorrelation
#'   panels (given as ORIGINAL parameter names). Default selects the
#'   population-level effect and between-study SD.
#' @param labels Optional named vector mapping original parameter names to
#'   display labels, e.g. \code{c(b_Intercept = "\u03bc", sd_Author__Intercept = "\u03c4")}.
#'   Overrides the automatic mu / tau defaults.
#' @param ndraws Draws used in the posterior predictive panel.
#' @param font Font family for panel text.
#' @param ... Unused.
#'
#' @return A patchwork object (the assembled diagnostic page).
#' @export
#'
#' @importFrom posterior as_draws_array summarise_draws ndraws variables
#'   subset_draws as_draws_matrix rhat ess_bulk ess_tail
diagnostics <- function(object, pars = NULL, labels = NULL, study_labels = NULL,
                        model_label = NULL, ndraws = 100, font = NULL, ...) {
  rlang::check_installed("bayesplot", reason = "to draw MCMC diagnostics.")
  UseMethod("diagnostics")
}

#' @export
diagnostics.default <- function(object, pars = NULL, labels = NULL, study_labels = NULL,
                                model_label = NULL, ndraws = 100, font = NULL, ...) {

  draws_array <- .diag_draws(object)
  summary_df  <- .diag_summary(object)

  # choose parameters on ORIGINAL names, then build the display-label map
  if (is.null(pars)) pars <- .diag_default_pars(draws_array)
  label_map <- .diag_param_labels(object, draws_array, labels, study_labels)

  # numeric health summary (uses original draws/summary)
  diag_stats <- .diag_extract(object, summary_df, draws_array)

  # relabel draws, summary and the chosen pars for display
  draws_disp <- .diag_relabel(draws_array, label_map)
  pars_disp  <- ifelse(pars %in% names(label_map), label_map[pars], pars)

  # drop the internal brms `Intercept` (centered-scale duplicate) from panels
  summary_df <- summary_df[summary_df$variable != "Intercept", , drop = FALSE]
  summary_df$variable <- ifelse(summary_df$variable %in% names(label_map),
                                label_map[summary_df$variable], summary_df$variable)

  # canonical display order for the R-hat / ESS panels:
  # mu, tau, nu, alpha, sigma, log prior, log posterior, then study names (A-Z)
  fixed_head <- c("\u03bc", "\u03c4", "\u03bd", "\u03b1", "\u03c3",
                  "log prior", "log posterior")
  present_head <- fixed_head[fixed_head %in% summary_df$variable]
  studies      <- sort(setdiff(summary_df$variable, present_head))
  param_order  <- c(present_head, studies)
  summary_df$variable <- factor(summary_df$variable, levels = rev(param_order))

  p_trace <- .diag_panel_trace(draws_disp, pars_disp, font)
  p_ppc   <- .diag_panel_ppc(object, ndraws, font)
  p_rhat  <- .diag_panel_rhat(summary_df, font)
  p_neff  <- .diag_panel_neff(summary_df, posterior::ndraws(draws_array), font)
  p_acf   <- .diag_panel_acf(draws_disp, pars_disp, font)
  p_table <- .diag_panel_table(diag_stats, font)

  diag_title <- if (is.null(model_label)) "Model Diagnostics"
  else paste0("Model Diagnostics: ", model_label)

  (p_trace | p_ppc) /
    (p_rhat  | p_neff) /
    (p_acf   | p_table) +
    patchwork::plot_annotation(
      title = diag_title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 16, family = font)))
}


# ---------------------------------------------------------------------------
# Backend accessors
# ---------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.diag_draws <- function(fit) UseMethod(".diag_draws")

#' @keywords internal
#' @noRd
.diag_draws.default     <- function(fit) posterior::as_draws_array(fit)

#' @keywords internal
#' @noRd
.diag_draws.CmdStanMCMC <- function(fit) posterior::as_draws_array(fit$draws())

#' @keywords internal
#' @noRd
.diag_summary <- function(fit) {
  posterior::summarise_draws(
    .diag_draws(fit),
    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk, ess_tail = posterior::ess_tail)
}

# population-level parameters first; exclude degenerate residual sigma (se()
# models), lp__, raw z's, per-study r_/theta[] and generated quantities.
#' @keywords internal
#' @noRd
.diag_default_pars <- function(draws_array) {
  all_pars  <- dimnames(draws_array)[[3]]
  preferred <- grep("^(b_Intercept|mu|pooled|sd_.*__Intercept|tau|nu|alpha)$",
                    all_pars, value = TRUE)
  preferred <- setdiff(preferred, "sigma")
  if (length(preferred) == 0) {
    preferred <- all_pars[!grepl(
      "^(lp__|lprior|sigma|z\\[|r_|theta\\[|shrinkage\\[|y_rep\\[|log_lik\\[)", all_pars)]
  }
  utils::head(preferred, 6)
}


# ---------------------------------------------------------------------------
# Display labels: pooled effect -> mu, between-study SD -> tau (+ nu / alpha)
# ---------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.diag_param_labels <- function(fit, draws_array, labels = NULL, study_labels = NULL) {
  v <- dimnames(draws_array)[[3]]

  mu_par  <- intersect(c("b_Intercept", "mu", "pooled"), v)[1]
  tau_par <- c(grep("^sd_.*__Intercept$", v, value = TRUE), intersect("tau", v))[1]
  nu_par  <- intersect("nu", v)[1]
  al_par  <- intersect("alpha", v)[1]

  m <- c(
    stats::setNames("\u03bc", mu_par),   # mu
    stats::setNames("\u03c4", tau_par),  # tau
    stats::setNames("\u03bd", nu_par),   # nu
    stats::setNames("\u03b1", al_par))   # alpha
  m <- m[!is.na(names(m))]

  # global / bookkeeping parameters (NOT the internal brms `Intercept`, which is
  # the centered-scale intercept duplicate of b_Intercept and should not be shown)
  bookkeeping <- c(lp__ = "log posterior", lprior = "log prior", sigma = "\u03c3")
  present_bk <- bookkeeping[names(bookkeeping) %in% v]
  m <- c(m, present_bk)

  # per-study random effects: r_<group>[<key>,Intercept] -> study display name
  re_terms <- grep("^r_.*\\[.*,Intercept\\]$", v, value = TRUE)
  if (length(re_terms) > 0) {
    keys <- sub("^r_[^\\[]*\\[(.*),Intercept\\]$", "\\1", re_terms)  # e.g. Ramanan_a
    disp <- if (!is.null(study_labels)) {
      # map de-duplicated key -> supplied display name; fall back to a tidy key
      out <- unname(study_labels[keys])
      out[is.na(out)] <- gsub("[._]", " ", keys[is.na(out)])
      out
    } else {
      gsub("[._]", " ", keys)   # tidy the raw key (Van.Zyl -> Van Zyl, Ramanan_a -> Ramanan a)
    }
    m <- c(m, stats::setNames(disp, re_terms))
  }

  if (!is.null(labels)) m[names(labels)] <- labels

  # guarantee unique display labels (bayesplot's mcmc_rhat/neff require it):
  # append a counter to any collisions, preserving the names (raw param -> label).
  if (any(duplicated(unname(m)))) {
    vals <- unname(m)
    dup  <- duplicated(vals) | duplicated(vals, fromLast = TRUE)
    vals[dup] <- make.unique(vals[dup], sep = " ")
    m <- stats::setNames(vals, names(m))
  }
  m
}

#' @keywords internal
#' @noRd
.diag_relabel <- function(draws_array, label_map) {
  if (length(label_map) == 0) return(draws_array)
  v   <- dimnames(draws_array)[[3]]
  hit <- v %in% names(label_map)
  v[hit] <- unname(label_map[v[hit]])
  dimnames(draws_array)[[3]] <- v
  draws_array
}


# ---------------------------------------------------------------------------
# NUTS diagnostics (divergences / treedepth / E-BFMI), per backend
# ---------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.diag_nuts <- function(fit) UseMethod(".diag_nuts")

#' @keywords internal
#' @noRd
.diag_nuts.CmdStanMCMC <- function(fit) {
  ds <- fit$diagnostic_summary()
  list(num_divergent = ds$num_divergent,
       num_max_treedepth = ds$num_max_treedepth,
       ebfmi = ds$ebfmi)
}

#' @keywords internal
#' @noRd
.diag_nuts.brmsfit <- function(fit) {
  np <- brms::nuts_params(fit)
  by_chain <- function(param, fun) {
    d <- np[np$Parameter == param, , drop = FALSE]
    if (nrow(d) == 0) return(NA_real_)
    tapply(d$Value, d$Chain, fun)
  }
  max_td <- tryCatch(fit$fit@stan_args[[1]]$control$max_treedepth,
                     error = function(e) NULL) %||% 10
  en <- np[np$Parameter == "energy__", , drop = FALSE]
  ebfmi <- if (nrow(en) == 0) NA_real_ else
    tapply(en$Value, en$Chain, function(e) sum(diff(e)^2) / sum((e - mean(e))^2))
  list(
    num_divergent     = as.integer(by_chain("divergent__", function(v) sum(v == 1))),
    num_max_treedepth = as.integer(by_chain("treedepth__", function(v) sum(v >= max_td))),
    ebfmi             = as.numeric(ebfmi))
}

#' @keywords internal
#' @noRd
.diag_nuts.stanfit <- function(fit) {
  rlang::check_installed("rstan")
  sp  <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  mtd <- tryCatch(fit@stan_args[[1]]$control$max_treedepth,
                  error = function(e) NULL) %||% 10
  list(
    num_divergent     = as.integer(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1))),
    num_max_treedepth = as.integer(vapply(sp, function(x) sum(x[, "treedepth__"] >= mtd), numeric(1))),
    ebfmi             = tryCatch(rstan::get_bfmi(fit), error = function(e) NA_real_))
}

#' @keywords internal
#' @noRd
.diag_nuts.default <- function(fit)
  list(num_divergent = NA_integer_, num_max_treedepth = NA_integer_, ebfmi = NA_real_)


# ---------------------------------------------------------------------------
# Numeric diagnostic summary
# ---------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.diag_extract <- function(object, summary_df, draws_array) {
  total_draws <- posterior::ndraws(draws_array)

  max_rhat    <- max(summary_df$rhat, na.rm = TRUE)
  n_high_rhat <- sum(summary_df$rhat > 1.01, na.rm = TRUE)

  ess_ratio <- summary_df$ess_bulk / total_draws
  min_ess   <- min(ess_ratio, na.rm = TRUE)
  n_low_ess <- sum(ess_ratio < 0.1, na.rm = TRUE)

  nuts   <- .diag_nuts(object)
  n_div  <- if (all(is.na(nuts$num_divergent)))     NA_integer_ else sum(nuts$num_divergent)
  max_td <- if (all(is.na(nuts$num_max_treedepth))) NA_integer_ else max(nuts$num_max_treedepth)
  minbf  <- if (all(is.na(nuts$ebfmi)))             NA_real_    else min(nuts$ebfmi, na.rm = TRUE)

  flags <- c(max_rhat <= 1.01,
             min_ess >= 0.1,
             if (is.na(n_div))  NA else n_div == 0,
             if (is.na(max_td)) NA else max_td == 0,
             if (is.na(minbf))  NA else minbf > 0.3)

  list(max_rhat = max_rhat, n_high_rhat = n_high_rhat,
       min_ess_ratio = min_ess, n_low_ess = n_low_ess,
       n_divergent = n_div, n_max_treedepth = max_td, min_ebfmi = minbf,
       all_ok = all(flags[!is.na(flags)]))
}


# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.diag_base_theme <- function(font = NULL) {
  ggplot2::theme_light() +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold", family = font),
      axis.title  = ggplot2::element_text(family = font),
      axis.text   = ggplot2::element_text(family = font),
      legend.text = ggplot2::element_text(family = font),
      panel.grid.minor = ggplot2::element_blank())
}

#' @keywords internal
#' @noRd
.diag_empty <- function(title, msg, font = NULL) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = msg, size = 3.2, family = font) +
    ggplot2::ggtitle(title) + ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", family = font))
}

#' @keywords internal
#' @noRd
.diag_panel_trace <- function(draws_array, pars, font = NULL) {
  avail <- intersect(pars, dimnames(draws_array)[[3]])
  if (length(avail) == 0) return(.diag_empty("Trace", "No parameters to display.", font))
  tryCatch(
    bayesplot::mcmc_trace(draws_array, pars = avail) +
      ggplot2::ggtitle("Trace") + .diag_base_theme(font),
    error = function(e) .diag_empty("Trace", conditionMessage(e), font))
}

#' @keywords internal
#' @noRd
.diag_panel_acf <- function(draws_array, pars, font = NULL) {
  avail <- intersect(pars, dimnames(draws_array)[[3]])
  if (length(avail) == 0) return(.diag_empty("Autocorrelation", "No parameters.", font))
  tryCatch(
    bayesplot::mcmc_acf(draws_array, pars = avail) +
      ggplot2::ggtitle("Autocorrelation") + .diag_base_theme(font),
    error = function(e) .diag_empty("Autocorrelation", conditionMessage(e), font))
}

#' @keywords internal
#' @noRd
.diag_panel_rhat <- function(summary_df, font = NULL) {
  d <- summary_df
  d <- d[!is.na(d$rhat), , drop = FALSE]
  d$band <- cut(d$rhat, breaks = c(-Inf, 1.01, 1.05, Inf),
                labels = c("\u2264 1.01", "\u2264 1.05", "> 1.05"))
  cols <- c("\u2264 1.01" = "#9ecae1", "\u2264 1.05" = "#4292c6", "> 1.05" = "#08306b")
  # fix x-range so tiny deviations near 1 don't get magnified (bayesplot behaviour)
  xmax <- max(1.05, max(d$rhat, na.rm = TRUE) * 1.001)
  tryCatch(
    ggplot2::ggplot(d, ggplot2::aes(x = rhat, y = variable, colour = band)) +
      ggplot2::geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.4) +
      ggplot2::geom_point(size = 2.4, show.legend = TRUE) +
      ggplot2::scale_colour_manual(values = cols, drop = FALSE,
                                   name = expression(hat(R)),
                                   na.translate = FALSE) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        override.aes = list(colour = unname(cols)))) +
      ggplot2::scale_x_continuous(limits = c(1, xmax)) +
      ggplot2::labs(x = expression(hat(R)), y = NULL, title = "R-hat") +
      .diag_base_theme(font),
    error = function(e) .diag_empty("R-hat", conditionMessage(e), font))
}

#' @keywords internal
#' @noRd
.diag_panel_neff <- function(summary_df, total_draws, font = NULL) {
  d <- summary_df
  d$ratio <- d$ess_bulk / total_draws
  d <- d[!is.na(d$ratio), , drop = FALSE]
  d$band <- cut(d$ratio, breaks = c(-Inf, 0.1, 0.5, Inf),
                labels = c("\u2264 0.1", "\u2264 0.5", "> 0.5"))
  cols <- c("\u2264 0.1" = "#08306b", "\u2264 0.5" = "#4292c6", "> 0.5" = "#9ecae1")
  xmax <- max(1, max(d$ratio, na.rm = TRUE) * 1.02)
  tryCatch(
    ggplot2::ggplot(d, ggplot2::aes(x = ratio, y = variable, colour = band)) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = ratio,
                                         y = variable, yend = variable),
                            colour = "grey85", linewidth = 0.4) +
      ggplot2::geom_point(size = 2.4, show.legend = TRUE) +
      ggplot2::scale_colour_manual(values = cols, drop = FALSE,
                                   name = expression(N[eff]/N),
                                   na.translate = FALSE) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        override.aes = list(colour = unname(cols)))) +
      ggplot2::scale_x_continuous(limits = c(0, xmax)) +
      ggplot2::labs(x = expression(N[eff]/N), y = NULL, title = "ESS ratio") +
      .diag_base_theme(font),
    error = function(e) .diag_empty("ESS ratio", conditionMessage(e), font))
}

#' @keywords internal
#' @noRd
.diag_panel_ppc <- function(object, ndraws = 100, font = NULL) {
  tryCatch({
    if (inherits(object, "brmsfit")) {
      return(brms::pp_check(object, ndraws = ndraws) +
               ggplot2::ggtitle("Posterior Predictive Check") + .diag_base_theme(font))
    }
    dr <- .diag_draws(object)
    yrep_vars <- grep("^y_rep(\\[|$)", posterior::variables(dr), value = TRUE)
    if (length(yrep_vars) == 0)
      return(.diag_empty("Posterior Predictive Check",
                         "No y_rep in draws.\nAdd a generated quantities block.", font))
    y_rep <- posterior::as_draws_matrix(posterior::subset_draws(dr, variable = yrep_vars))
    if (nrow(y_rep) > ndraws) y_rep <- y_rep[sample.int(nrow(y_rep), ndraws), , drop = FALSE]
    y_obs <- tryCatch(object$data$y %||% object$data$yi, error = function(e) NULL)
    if (is.null(y_obs))
      return(.diag_empty("Posterior Predictive Check",
                         "Observed outcome not stored in fit.", font))
    bayesplot::ppc_dens_overlay(y = y_obs, yrep = y_rep) +
      ggplot2::ggtitle("Posterior Predictive Check") + .diag_base_theme(font)
  }, error = function(e) .diag_empty("Posterior Predictive Check", conditionMessage(e), font))
}

# ---- pass/fail table: coloured fill + plain words (device-independent) -------
#' @keywords internal
#' @noRd
.diag_panel_table <- function(d, font = NULL) {

  num <- function(v, digits) if (is.na(v)) "\u2014" else formatC(v, format = "f", digits = digits)
  int <- function(v)        if (is.na(v)) "\u2014" else as.character(v)

  rows <- tibble::tibble(
    y     = 6:1,
    Check = c("Max R-hat", "Min ESS ratio", "Divergences",
              "Max-treedepth hits", "Min E-BFMI", "Overall"),
    val   = c(num(d$max_rhat, 3),
              num(d$min_ess_ratio, 2),
              int(d$n_divergent),
              int(d$n_max_treedepth),
              num(d$min_ebfmi, 2),
              if (isTRUE(d$all_ok)) "PASS" else "REVIEW"),
    ok    = c(d$max_rhat <= 1.01,
              d$min_ess_ratio >= 0.1,
              if (is.na(d$n_divergent)) NA else d$n_divergent == 0,
              if (is.na(d$n_max_treedepth)) NA else d$n_max_treedepth == 0,
              if (is.na(d$min_ebfmi)) NA else d$min_ebfmi > 0.3,
              d$all_ok)) |>
    dplyr::mutate(
      fill = dplyr::case_when(is.na(ok) ~ "grey90",
                              ok        ~ "#d9ead3",   # green = pass
                              TRUE      ~ "#f4cccc"))   # red   = fail

  ggplot2::ggplot(rows, ggplot2::aes(y = y)) +
    ggplot2::geom_tile(ggplot2::aes(x = 1.55), width = 0.85, height = 0.82,
                       fill = rows$fill, colour = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(x = 0, label = Check),
                       hjust = 0, family = font, size = 3.6) +
    ggplot2::geom_text(ggplot2::aes(x = 1.55, label = val),
                       family = font, fontface = "bold", size = 3.6) +
    ggplot2::scale_x_continuous(limits = c(-0.05, 2.1)) +
    ggplot2::scale_y_continuous(limits = c(0.4, 6.6)) +
    ggplot2::ggtitle("MCMC Health") +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", family = font),
      plot.margin = ggplot2::margin(5, 5, 5, 5))
}

# null-coalesce (define only if not already available)
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a


# ---------------------------------------------------------------------------
# Helper: build a study_labels map (de-duplicated key -> "Author (year)")
# for the per-study random-effect labels in the R-hat / ESS panels.
#
# `data` is the analysis dataset with the (make_authors_unique) Author key,
# plus the original label and year. Pass the result as
#   diagnostics(fit, study_labels = make_study_labels(dat_primary))
# ---------------------------------------------------------------------------
#' Build display labels for per-study random effects
#'
#' @param data Analysis dataset containing the (unique) \code{Author} key and,
#'   ideally, \code{Author_original} and \code{year}.
#' @return A named character vector: names = Author keys (as they appear in
#'   \code{r_Author[<key>,Intercept]}); values = "Author (year)".
#' @export
make_study_labels <- function(data) {
  key <- as.character(data$Author)
  base <- if ("Author_original" %in% names(data))
    as.character(data$Author_original) else gsub("[._]", " ", key)
  yr  <- if ("year" %in% names(data)) data$year else NA
  disp <- ifelse(is.na(yr), base, sprintf("%s (%s)", base, yr))
  stats::setNames(disp, key)
}
