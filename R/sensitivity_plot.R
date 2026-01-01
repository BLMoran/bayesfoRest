#' Generate Sensitivity Analysis Plot for Bayesian Meta-Analysis
#'
#' @description
#' Creates a comprehensive sensitivity analysis visualization showing how meta-analytic 
#' estimates vary across different prior specifications and analysis approaches. The plot
#' combines forest-plot style density distributions with tabulated results.
#'
#' @param model A fitted brmsfit object from a Bayesian meta-analysis
#' @param data A data frame containing the study data used to fit the model
#' @param priors A named list of brms prior specifications (created with \code{brms::prior()})
#'   with names typically including "vague", "weakreg", and "informative"
#' @param measure Character string specifying the effect measure. One of "OR" (odds ratio),
#'   "RR" (risk ratio), "HR" (hazard ratio), "IRR" (incidence rate ratio), 
#'   "MD" (mean difference), or "SMD" (standardized mean difference)
#' @param study_var Optional. Name of the study identifier variable (unquoted). Required
#'   if \code{incl_pet_peese = TRUE} or \code{incl_mixture = TRUE}
#' @param rob_var Optional. Name of the risk of bias variable (unquoted). Should contain
#'   values like "Low", "Some concerns", "High"
#' @param exclude_high_rob Logical. If TRUE, performs sensitivity analysis excluding
#'   studies with high risk of bias. Default is FALSE
#' @param incl_pet_peese Logical. If TRUE, includes PET-PEESE bias adjustment analysis.
#'   Default is FALSE
#' @param pet_peese_direction Character string. Direction of expected bias for PET-PEESE.
#'   Either "negative" or "positive". Default is "negative"
#' @param pet_peese_threshold Numeric. Threshold for choosing between PET and PEESE 
#'   (based on proportion of studies showing effect). Default is 0.10
#' @param incl_mixture Logical. If TRUE, includes robust mixture model analysis.
#'   Default is FALSE
#' @param incl_bma Logical. If TRUE, includes Bayesian model averaging results.
#'   Default is FALSE
#' @param model_bma Named list of RoBMA model fits (from \code{RoBMA_brms()}) to include
#'   in model averaging. Names should be "vague", "weakreg", or "informative"
#' @param add_probs Logical. If TRUE, adds probability columns to the results table.
#'   Default is FALSE
#' @param null_value Numeric. The null value for the effect measure. If NULL, uses
#'   default based on measure (1 for ratio measures, 0 for difference measures)
#' @param null_range Numeric vector of length 2 or single numeric value. Defines the
#'   range of practical equivalence around the null value. If single value, creates
#'   symmetric range around null_value
#' @param add_null_range Logical. If TRUE, displays the null range region on the plot.
#'   Default is FALSE
#' @param color_null_range Character string. Color for the null range region.
#'   Default is "#77bb41" (green)
#' @param label_control Character string. Label for the control/reference group.
#'   Default is "Control"
#' @param label_intervention Character string. Label for the intervention group.
#'   Default is "Intervention"
#' @param title Optional character string. Main title for the plot
#' @param subtitle Optional character string. Subtitle for the plot
#' @param title_align Character string. Alignment for title and subtitle.
#'   One of "left", "center", or "right". Default is "left"
#' @param xlim Optional numeric vector of length 2. X-axis limits for the density plot
#' @param x_breaks Optional numeric vector. Custom x-axis break points
#' @param color_palette Optional color palette for the plot
#' @param color_overall_posterior Character string. Fill color for posterior densities.
#'   Default is "dodgerblue"
#' @param color_overall_posterior_outline Character string. Outline color for posterior
#'   densities. Default is "blue"
#' @param split_color_by_null Logical. If TRUE, colors densities based on which side
#'   of null they favor. Default is FALSE
#' @param color_favours_control Character string. Color when effect favors control.
#'   Default is "firebrick"
#' @param color_favours_intervention Character string. Color when effect favors 
#'   intervention. Default is "dodgerblue"
#' @param plot_width Numeric. Width ratio for the density plot section. Default is 4
#' @param font Optional character string. Font family to use for the plot
#'
#' @return A patchwork object combining the sensitivity analysis tables and density plots
#'
#' @details
#' This function performs comprehensive sensitivity analyses including:
#' \itemize{
#'   \item Prior sensitivity: Re-fits the model with different prior specifications
#'   \item Risk of bias sensitivity: Excludes high risk of bias studies
#'   \item PET-PEESE: Publication bias adjustment using precision-based methods
#'   \item Mixture models: Robust analysis using mixture distributions
#'   \item Bayesian model averaging: Combines results across model specifications
#' }
#'
#' The output combines three elements using patchwork:
#' \itemize{
#'   \item Left table: Prior specifications
#'   \item Center plot: Density distributions of posterior estimates
#'   \item Right table: Numerical summaries and probabilities
#' }
#'
#' @examples
#' \dontrun{
#' # Basic sensitivity plot
#' sensitivity_plot(
#'   model = my_brms_model,
#'   data = my_data,
#'   priors = list(
#'     vague = prior("normal(0, 10)", class = "Intercept"),
#'     weakreg = prior("normal(0, 1)", class = "Intercept"),
#'     informative = prior("normal(0.5, 0.5)", class = "Intercept")
#'   ),
#'   measure = "SMD"
#' )
#'
#' # Advanced sensitivity with all options
#' sensitivity_plot(
#'   model = my_model,
#'   data = my_data,
#'   priors = my_priors,
#'   measure = "OR",
#'   study_var = study_id,
#'   rob_var = rob,
#'   exclude_high_rob = TRUE,
#'   incl_pet_peese = TRUE,
#'   incl_mixture = TRUE,
#'   incl_bma = TRUE,
#'   model_bma = my_robma_models,
#'   add_probs = TRUE,
#'   null_range = c(0.9, 1.1),
#'   add_null_range = TRUE
#' )
#' }
#'
#' @export
#' 
sensitivity_plot <- function(model,
                             data,
                             priors,
                             measure,
                             study_var = NULL,
                             rob_var = NULL,
                             exclude_high_rob = FALSE,
                             incl_pet_peese = FALSE,
                             pet_peese_direction = "negative",
                             pet_peese_threshold = 0.10,
                             incl_mixture = FALSE,
                             incl_bma = FALSE,
                             model_bma = NULL,
                             add_probs = FALSE,
                             null_value = NULL,
                             null_range = NULL,
                             add_null_range = FALSE,
                             color_null_range = "#77bb41",
                             label_control = "Control",
                             label_intervention = "Intervention",
                             title = NULL,
                             subtitle = NULL,
                             title_align = "left",
                             xlim = NULL,
                             x_breaks = NULL,
                             color_palette = NULL,
                             color_overall_posterior = "dodgerblue",
                             color_overall_posterior_outline = "blue",
                             split_color_by_null = FALSE,
                             color_favours_control = "firebrick",
                             color_favours_intervention = "dodgerblue",
                             plot_width = 4,
                             font = NULL) {
  
  # ---------------------------
  # 1. Validation
  # ---------------------------
  validate_inputs_sens_plot(
    model   = model,
    data    = data,
    measure = measure,
    priors = priors
  )
  
  if (incl_bma && is.null(model_bma)) {
    stop("`model_bma` must be supplied when `incl_bma = TRUE`.", call. = FALSE)
  }
  
  # Capture tidy eval variables
  rob_var <- rlang::enquo(rob_var)
  study_var <- rlang::enquo(study_var)
  
  # Validate study_var if PET-PEESE is included
  if (isTRUE(incl_pet_peese) && rlang::quo_is_null(study_var)) {
    stop("`study_var` must be supplied when `incl_pet_peese = TRUE`.", call. = FALSE)
  }
  
  props <- get_measure_properties(measure)
  null_value <- null_value %||% props$null_value
  
  # ---------------------------
  # 2. Null range handling
  # ---------------------------
  if (is.null(null_range) && isTRUE(add_null_range)) {
    null_range <- switch(
      measure,
      OR  = c(0.9, 1.1),
      RR  = c(0.9, 1.1),
      HR  = c(0.9, 1.1),
      IRR = c(0.9, 1.1),
      SMD = c(-0.1, 0.1),
      stop("For MD, `null_range` must be supplied.", call. = FALSE)
    )
  } else if (!is.null(null_range)) {
    null_range <- if (length(null_range) == 1) {
      c(null_value - null_range, null_value + null_range)
    } else {
      sort(null_range)
    }
  }
  
  # ---------------------------
  # 3. Prior labels
  # ---------------------------
  prior_labels <- c(
    vague        = "Vague",
    weakreg     = "Weakly Regularising",
    informative = "Informative"
  )
  
  # ---------------------------
  # 4. Build draws (ALL SECTIONS)
  # ---------------------------
  sens_specs <- list(
    list(
      id = "all",
      include = TRUE,
      model = model,
      data = data,
      section_label = "All studies"
    ),
    list(
      id = "excl_high_rob",
      include = isTRUE(exclude_high_rob),
      model = model,
      data = dplyr::filter(data, !!rob_var != "High" | is.na(!!rob_var)),
      section_label = "Excluding High RoB"
    ),
    list(
      id = "pet_peese",
      include = isTRUE(incl_pet_peese),
      model = model,
      data = data,
      section_label = "PET / PEESE",
      study_var = study_var,
      pet_peese_direction = pet_peese_direction,
      pet_peese_threshold = pet_peese_threshold
    ),
    list(
      id = "mixture",
      include = isTRUE(incl_mixture),
      model = model,
      data = data,
      section_label = "Robust mixture model",
      study_var = study_var
    )
  )
  
  draws_ma <- sens_specs |>
    purrr::keep(~ .x$include) |>
    purrr::map_dfr(run_ma_sensitivity,
                   priors = priors,
                   measure = measure,
                   prior_labels = prior_labels)
  
  if (isTRUE(incl_bma)) {
    draws_bma <- purrr::imap_dfr(
      model_bma,
      ~ robma_to_sensitivity_draws(
        robma_fit   = .x,
        prior       = .y,
        prior_label = c(
          vague = "Vague",
          weakreg = "Weakly Regularising",
          informative = "Informative"
        )[[.y]],
        section_label = "Model-averaged"))
  } else {
    draws_bma <- NULL
  }
  
  draws <- dplyr::bind_rows(draws_ma, draws_bma)
  
  # ---------------------------
  # 5. Prior table
  # ---------------------------
  
  prior_table_ma <- extract_mu_tau_priors(
    priors = priors,
    model  = model)
  
  # Initialize empty tibble
  prior_table_bma <- tibble::tibble(
    prior = character(),
    prior_label = character(),
    mu_prior_unicode = character(),
    tau_prior_unicode = character()
  )
  
  if (incl_bma) {
    prior_table_bma <- purrr::imap_dfr(
      model_bma,
      ~ {
        priors <- extract_priors_from_robma_fit(.x)
        
        tibble::tibble(
          prior = .y,
          prior_label = c(
            vague = "Vague",
            weakreg = "Weakly Regularising",
            informative = "Informative"
          )[[.y]],
          mu_prior_unicode  = priors$mu_prior_unicode,
          tau_prior_unicode = priors$tau_prior_unicode
        )
      }
    )
    
    prior_table <- dplyr::bind_rows(
      extract_mu_tau_priors(priors = priors, model = model),
      prior_table_bma
    )
  } else {
    prior_table <- extract_mu_tau_priors(priors = priors, model = model)
  }
  
  # ---------------------------
  # 6. Posterior summary + join
  # ---------------------------
  
  summary <- draws |>
    summarise_sensitivity_posteriors(
      null_value = null_value,
      null_range = null_range
    ) |>
    dplyr::mutate(
      mu_prior_unicode = dplyr::if_else(
        section_label == "Model-averaged",
        prior_table_bma$mu_prior_unicode[
          match(paste(prior, prior_label),
                paste(prior_table_bma$prior, prior_table_bma$prior_label))
        ],
        prior_table_ma$mu_prior_unicode[
          match(paste(prior, prior_label),
                paste(prior_table_ma$prior, prior_table_ma$prior_label))
        ]
      ),
      tau_prior_unicode = dplyr::if_else(
        section_label == "Model-averaged",
        prior_table_bma$tau_prior_unicode[
          match(paste(prior, prior_label),
                paste(prior_table_bma$prior, prior_table_bma$prior_label))
        ],
        prior_table_ma$tau_prior_unicode[
          match(paste(prior, prior_label),
                paste(prior_table_ma$prior, prior_table_ma$prior_label))
        ]
      )
    )
  
  # ---------------------------
  # 7. Tables
  # ---------------------------
  
  table_left <- sensitivity_table_left(summary, font = font)
  table_right <- sensitivity_table_right(
    summary,
    measure   = measure,
    add_probs = add_probs,
    font      = font
  )
  
  # ---------------------------
  # 8. Density plot
  # ---------------------------
  
  sens_plot <- sensitivity_density_plot_fn(
    df = draws,
    measure = measure,
    split_color_by_null = split_color_by_null,
    color_overall_posterior = color_overall_posterior,
    color_overall_posterior_outline = color_overall_posterior_outline,
    color_favours_control = color_favours_control,
    color_favours_intervention = color_favours_intervention,
    label_control = label_control,
    label_intervention = label_intervention,
    xlim = xlim,
    x_breaks = x_breaks,
    null_value = null_value,
    null_range = null_range,
    color_null_range = color_null_range,
    add_null_range = add_null_range,
    font = font
  )
  
  # ---------------------------
  # 9. Patchwork assembly
  # ---------------------------
  
  patchwork_fn(
    table.left        = table_left,
    study.density.plot = sens_plot,
    table.right       = table_right,
    plot_width        = plot_width,
    title             = title,
    subtitle          = subtitle,
    title_align       = title_align,
    add_rob_legend    = FALSE,
    rob_tool          = "rob2",
    font              = font
  )
}
