#' Validate inputs for bayes_forest (updated for multi-model support)
#'
#' @noRd
validate_inputs <- function(
    model,
    data,
    measure,
    studyvar,
    year,
    subgroup,
    subgroup_var,
    sort_studies_by,
    shrinkage_output
) {
  # Validate model using new system
  validate_model(model)
  
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.")
  }
  
  rlang::arg_match(measure, c("OR", "HR", "RR", "IRR", "MD", "SMD"))
  rlang::arg_match(sort_studies_by, c("author", "year", "effect"))
  rlang::arg_match(shrinkage_output, c("density", "pointinterval"))
  
  if (rlang::quo_is_null(rlang::enquo(studyvar))) {
    rlang::abort("`studyvar` must be provided (bare column name).")
  }
  
  studyvar_name <- rlang::as_name(rlang::ensym(studyvar))
  if (!studyvar_name %in% names(data)) {
    rlang::abort(paste0("Can't find column `", studyvar_name, "` in `data`."))
  }
  
  if (!rlang::quo_is_null(rlang::enquo(year))) {
    year_name <- rlang::as_name(rlang::ensym(year))
    if (!year_name %in% names(data)) {
      rlang::abort(paste0("Can't find column `", year_name, "` in `data`."))
    }
  }
  
  # Validate subgroup variable if subgroup = TRUE
  if (isTRUE(subgroup)) {
    if (rlang::quo_is_null(rlang::enquo(subgroup_var))) {
      if (!"Subgroup" %in% names(data)) {
        rlang::abort("`subgroup = TRUE` but no `subgroup_var` specified and no 'Subgroup' column found in data.")
      }
    } else {
      subgroup_var_name <- rlang::as_name(rlang::ensym(subgroup_var))
      if (!subgroup_var_name %in% names(data)) {
        rlang::abort(paste0("Can't find subgroup column `", subgroup_var_name, "` in `data`."))
      }
    }
  }
  
  invisible(TRUE)
}

#' Validate Model for bayesfoRest
#'
#' @description
#' Checks that a model object is supported and contains the necessary
#' components for forest plot generation.
#'
#' @param model A fitted model object
#'
#' @return Invisible TRUE if valid, otherwise throws an error
#'
#' @export
validate_model <- function(model) {
  model_type <- detect_model_type(model)
  
  if (model_type == "unknown") {
    rlang::abort(
      c(
        paste0("Unsupported model type: '", class(model)[1], "'"),
        "i" = "Supported types: brmsfit, meta_stan, bayesmeta, stanreg, stanfit, CmdStanMCMC",
        "i" = "Use `detect_model_type(model)` to check model type."
      )
    )
  }
  
  # Type-specific validation
  if (model_type == "metastan") {
    if (is.null(model$fit)) {
      rlang::abort("MetaStan model does not contain fitted Stan object.")
    }
  }
  
  if (model_type == "rstanarm") {
    if (is.null(model$stanfit)) {
      rlang::abort("rstanarm model does not contain fitted Stan object.")
    }
  }
  
  if (model_type == "cmdstanr") {
    if (!model$runset$completed()) {
      rlang::abort("CmdStanR model has not completed sampling.")
    }
  }
  
  invisible(TRUE)
}


#' Validate Inputs for Sensitivity Plot
#'
#' @description
#' Internal validation function that checks the validity of core inputs to the
#' sensitivity plot function.
#'
#' @param model A brmsfit object to validate
#' @param data A data frame to validate
#' @param measure Character string specifying the effect measure
#'
#' @return Invisible TRUE if validation passes, otherwise throws an error
#'
#' @keywords internal
#' @noRd
validate_inputs_sens_plot <- function(
    model,
    data,
    priors,
    measure
) {
  
  if (!inherits(model, "brmsfit")) {
    rlang::abort("`model` must be a brmsfit object.")
  }
  
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.")
  }
  
  rlang::arg_match(measure, c("OR", "HR", "RR", "IRR", "MD", "SMD"))
  
  validate_priors(priors)
  
  invisible(TRUE)
}

#' Validate Priors Object
#'
#' @description
#' Checks that priors is either a single brms prior or a list of brms priors.
#'
#' @param priors Object to validate as brms prior(s)
#'
#' @return Invisible TRUE if valid, otherwise throws an error
#'
#' @keywords internal
#' @noRd
validate_priors <- function(priors) {
  
  if (is.null(priors)) {
    return(invisible(TRUE))
  }
  
  if (inherits(priors, "brmsprior")) {
    return(invisible(TRUE))
  }
  
  if (is.list(priors) && all(vapply(priors, inherits, logical(1), "brmsprior"))) {
    return(invisible(TRUE))
  }
  
  rlang::abort(
    "`priors` must be a brms prior or a list of brms priors created with brms::prior()."
  )
}

#' Ensure yi and vi Columns Exist in Data
#'
#' @description
#' For models that don't use metafor's escalc() (like MetaStan with raw counts),
#' this function calculates yi (log OR) and vi (variance) from the raw data
#' if they don't already exist.
#'
#' @param data Data frame with study information
#' @param measure Effect measure type ("OR", "RR", etc.)
#'
#' @return Data frame with yi and vi columns added if missing
#'
#' @keywords internal
#' @noRd
ensure_yi_vi <- function(data, measure = "OR") {
  
  # If yi and vi already exist, return as-is
  if (all(c("yi", "vi") %in% names(data))) {
    return(data)
  }
  
  # Calculate based on measure type and available columns
  if (measure %in% c("OR", "RR")) {
    # Need event counts and sample sizes
    # Check for various column naming conventions
    has_events <- all(c("Event_Intervention", "Event_Control") %in% names(data)) ||
      all(c("r1", "r2") %in% names(data)) ||
      all(c("rt", "rc") %in% names(data))
    
    has_sizes <- all(c("N_Intervention", "N_Control") %in% names(data)) ||
      all(c("n1", "n2") %in% names(data)) ||
      all(c("nt", "nc") %in% names(data))
    
    if (has_events && has_sizes) {
      # Get column names
      if (all(c("Event_Intervention", "Event_Control") %in% names(data))) {
        r1_col <- "Event_Intervention"; r2_col <- "Event_Control"
        n1_col <- "N_Intervention"; n2_col <- "N_Control"
      } else if (all(c("r1", "r2") %in% names(data))) {
        r1_col <- "r1"; r2_col <- "r2"
        n1_col <- "n1"; n2_col <- "n2"
      } else {
        r1_col <- "rt"; r2_col <- "rc"
        n1_col <- "nt"; n2_col <- "nc"
      }
      
      data <- data |>
        dplyr::mutate(
          # Apply 0.5 continuity correction if needed
          .r1_adj = dplyr::if_else(
            .data[[r1_col]] == 0 | .data[[r2_col]] == 0 |
              (.data[[n1_col]] - .data[[r1_col]]) == 0 |
              (.data[[n2_col]] - .data[[r2_col]]) == 0,
            .data[[r1_col]] + 0.5, as.numeric(.data[[r1_col]])
          ),
          .r2_adj = dplyr::if_else(
            .data[[r1_col]] == 0 | .data[[r2_col]] == 0 |
              (.data[[n1_col]] - .data[[r1_col]]) == 0 |
              (.data[[n2_col]] - .data[[r2_col]]) == 0,
            .data[[r2_col]] + 0.5, as.numeric(.data[[r2_col]])
          ),
          .n1_adj = dplyr::if_else(
            .data[[r1_col]] == 0 | .data[[r2_col]] == 0 |
              (.data[[n1_col]] - .data[[r1_col]]) == 0 |
              (.data[[n2_col]] - .data[[r2_col]]) == 0,
            .data[[n1_col]] + 1, as.numeric(.data[[n1_col]])
          ),
          .n2_adj = dplyr::if_else(
            .data[[r1_col]] == 0 | .data[[r2_col]] == 0 |
              (.data[[n1_col]] - .data[[r1_col]]) == 0 |
              (.data[[n2_col]] - .data[[r2_col]]) == 0,
            .data[[n2_col]] + 1, as.numeric(.data[[n2_col]])
          )
        )
      
      if (measure == "OR") {
        # Log odds ratio
        data <- data |>
          dplyr::mutate(
            yi = log((.r1_adj / (.n1_adj - .r1_adj)) / (.r2_adj / (.n2_adj - .r2_adj))),
            vi = 1/.r1_adj + 1/(.n1_adj - .r1_adj) + 1/.r2_adj + 1/(.n2_adj - .r2_adj)
          )
      } else if (measure == "RR") {
        # Log risk ratio
        data <- data |>
          dplyr::mutate(
            yi = log((.r1_adj / .n1_adj) / (.r2_adj / .n2_adj)),
            vi = 1/.r1_adj - 1/.n1_adj + 1/.r2_adj - 1/.n2_adj
          )
      }
      
      # Remove temporary columns
      data <- data |>
        dplyr::select(-dplyr::starts_with(".r"), -dplyr::starts_with(".n"))
      
    } else {
      rlang::warn(
        c(
          "Cannot calculate yi/vi: missing required columns for binary outcome.",
          "i" = "Need event counts (Event_Intervention/Event_Control or r1/r2) and",
          "i" = "sample sizes (N_Intervention/N_Control or n1/n2).",
          "i" = "Shrinkage plots will not display correctly."
        )
      )
      data$yi <- NA_real_
      data$vi <- NA_real_
    }
    
  } else if (measure %in% c("MD", "SMD")) {
    # Need means, SDs, and sample sizes
    has_means <- all(c("Mean_Intervention", "Mean_Control") %in% names(data))
    has_sds <- all(c("SD_Intervention", "SD_Control") %in% names(data))
    has_sizes <- all(c("N_Intervention", "N_Control") %in% names(data))
    
    if (has_means && has_sds && has_sizes) {
      if (measure == "MD") {
        data <- data |>
          dplyr::mutate(
            yi = Mean_Intervention - Mean_Control,
            vi = (SD_Intervention^2 / N_Intervention) + (SD_Control^2 / N_Control)
          )
      } else if (measure == "SMD") {
        # Hedges' g approximation
        data <- data |>
          dplyr::mutate(
            .pooled_sd = sqrt(((N_Intervention - 1) * SD_Intervention^2 + 
                                 (N_Control - 1) * SD_Control^2) / 
                                (N_Intervention + N_Control - 2)),
            .j = 1 - 3 / (4 * (N_Intervention + N_Control - 2) - 1),
            yi = .j * (Mean_Intervention - Mean_Control) / .pooled_sd,
            vi = (N_Intervention + N_Control) / (N_Intervention * N_Control) +
              yi^2 / (2 * (N_Intervention + N_Control))
          ) |>
          dplyr::select(-dplyr::starts_with("."))
      }
    } else {
      rlang::warn(
        c(
          "Cannot calculate yi/vi: missing required columns for continuous outcome.",
          "i" = "Need Mean_Intervention, Mean_Control, SD_Intervention, SD_Control,",
          "i" = "N_Intervention, and N_Control.",
          "i" = "Shrinkage plots will not display correctly."
        )
      )
      data$yi <- NA_real_
      data$vi <- NA_real_
    }
    
  } else {
    # For other measures, set to NA
    rlang::warn(
      c(
        paste0("Cannot calculate yi/vi for measure '", measure, "'."),
        "i" = "Please provide yi and vi columns in your data.",
        "i" = "Shrinkage plots will not display correctly."
      )
    )
    data$yi <- NA_real_
    data$vi <- NA_real_
  }
  
  return(data)
}