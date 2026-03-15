#' Internal function to disambiguate duplicate Author names
#'
#' When multiple studies share the same Author name, this function appends a
#' letter suffix to make each row unique. Two layers of disambiguation are
#' applied:
#'
#' 1. **Display-level (Year suffix):** When the same Author appears with the
#'    same Year (e.g. multi-arm trials), a letter is appended to the Year
#'    column (e.g. "1999" -> "1999a", "1999b"). This ensures the tables show
#'    "Smith (1999a)" and "Smith (1999b)".
#'
#' 2. **Internal-level (Author suffix):** When the same Author name appears
#'    more than once (regardless of Year), the Author column itself is made
#'    unique by appending a letter (e.g. "Smith" -> "Smith_a", "Smith_b").
#'    This is needed so that brms sees unique grouping levels and
#'    ggplot2 maps each study to its own row on the y-axis.
#'
#' The original Author name is preserved in a new column `Author_original` so
#' it can be restored for display purposes later.
#'
#' @param data A data frame containing at least `Author` and `Year` columns.
#' @return The data frame with unique `Author` values, disambiguated `Year`
#'   values where needed, and a new `Author_original` column containing the
#'   original (possibly duplicated) names.
#' @noRd
make_authors_unique <- function(data) {
  # Preserve the original name for display in tables
  data$Author_original <- data$Author
  
  # --- Step 1: Disambiguate Year for same Author + same Year combos ---
  # This ensures table display shows e.g. "Smith (1999a)", "Smith (1999b)"
  # Convert Year to character so we can append letter suffixes
  data$Year <- as.character(data$Year)
  
  author_year_combos <- paste(data$Author, data$Year, sep = "|||")
  dup_combos <- author_year_combos[duplicated(author_year_combos)]
  
  if (length(dup_combos) > 0) {
    for (combo in unique(dup_combos)) {
      idx <- which(author_year_combos == combo)
      suffixes <- letters[seq_along(idx)]
      # Extract the original year from the combo
      original_year <- strsplit(combo, "|||", fixed = TRUE)[[1]][2]
      data$Year[idx] <- paste0(original_year, suffixes)
    }
  }
  
  # --- Step 2: Disambiguate Author names that appear more than once ---
  # This is needed regardless of Year, because brms and ggplot2 use Author
  # as a grouping / factor variable and need unique levels.
  dup_authors <- data$Author[duplicated(data$Author)]
  
  if (length(dup_authors) > 0) {
    for (author_name in unique(dup_authors)) {
      idx <- which(data$Author == author_name)
      suffixes <- letters[seq_along(idx)]
      data$Author[idx] <- paste0(author_name, "_", suffixes)
    }
  }
  
  return(data)
}

#' Internal function to sort studies
#'
#' @noRd
sort_studies_fn <- function(.data, sort_studies_by = NULL) {
  sort_by <- if (is.null(sort_studies_by)) "author" else sort_studies_by
  has_subgroup <- "Subgroup" %in% names(.data) &&
    !all(is.na(.data$Subgroup)) &&
    dplyr::n_distinct(.data$Subgroup) > 1
  
  # Separate rows (Pooled Effect & Prediction) from actual studies.
  pooled <- .data |> dplyr::filter(Author == "Pooled Effect")
  prediction <- .data |> dplyr::filter(Author == "Prediction")
  studies <- .data |> dplyr::filter(!Author %in% c("Pooled Effect", "Prediction"))
  
  # Arrange studies
  studies <- switch(
    sort_by,
    "author" = studies |> dplyr::arrange(dplyr::across(any_of("Subgroup")), Author),
    "year" = studies |> dplyr::arrange(dplyr::across(any_of("Subgroup")), Year, Author),
    "effect" = studies |> dplyr::arrange(dplyr::across(any_of("Subgroup")), yi, Author),
    stop("Invalid value for `sort_studies_by`. Must be one of 'author', 'year', or 'effect'.")
  )
  
  # Recombine: studies, then Pooled Effect, then Prediction
  full <- dplyr::bind_rows(studies, pooled, prediction)
  
  full <- full |>
    dplyr::mutate(
      Author = factor(Author, levels = unique(Author))
    )
  
  return(full)
}

#' Get subgroup order based on sorting preference
#'
#' @noRd
get_subgroup_order <- function(data, sort_subgroup_by) {
  # Early return if no subgroup column
  if (!"Subgroup" %in% names(data)) {
    return(NULL)
  }
  
  # Get unique non-NA subgroups
  unique_subgroups <- unique(data$Subgroup[!is.na(data$Subgroup)])
  
  # Handle different sorting options
  if (is.character(sort_subgroup_by)) {
    if (length(sort_subgroup_by) == 1) {
      # Single string: either "alphabetical" or "effect"
      subgroup_order <- switch(
        sort_subgroup_by,
        "alphabetical" = sort(unique_subgroups),
        "effect" = {
          data |>
            dplyr::group_by(Subgroup) |>
            dplyr::summarise(mean_effect = mean(yi, na.rm = TRUE), .groups = "drop") |>
            dplyr::arrange(mean_effect) |>
            dplyr::pull(Subgroup)
        },
        rlang::abort(
          paste0("Invalid sort_subgroup_by value: '", sort_subgroup_by,
                 "'. Must be 'alphabetical', 'effect', or a character vector of subgroup names.")
        )
      )
    } else {
      # Multiple strings: custom order provided by user
      # Validate that provided subgroups exist in data
      missing_groups <- setdiff(sort_subgroup_by, unique_subgroups)
      if (length(missing_groups) > 0) {
        rlang::warn(
          paste0("These subgroups in sort_subgroup_by were not found in data: ",
                 paste(missing_groups, collapse = ", "))
        )
      }
      subgroup_order <- sort_subgroup_by
    }
  } else {
    rlang::abort(
      "sort_subgroup_by must be a character vector. ",
      "Use 'alphabetical', 'effect', or provide a vector of subgroup names."
    )
  }
  
  # Always append "Overall" at the end
  c(subgroup_order, "Overall")
}


# =============================================================================
# MODEL-SPECIFIC UPDATE FUNCTIONS
# =============================================================================
#' Update Model with New Data (Generic)
#'
#' @description
#' S3 generic for updating a model with new/modified data. This is needed
#' when Author names are disambiguated or for subgroup analysis.
#'
#' @param model A fitted model object
#' @param data New data frame
#' @param ... Additional arguments
#'
#' @return Updated model object, or original if update not supported
#'
#' @keywords internal
#' @noRd
update_model <- function(model, data, ...) {
  UseMethod("update_model")
}

#' @export
update_model.default <- function(model, data, ...) {
  # For models that don't support updating, return original
  # The adapter will handle the mapping
  model
}

#' @export
update_model.brmsfit <- function(model, data, studyvar = Author, ...) {
  studyvar_name <- rlang::as_name(rlang::ensym(studyvar))
  
  # Extract response variable from formula
  full_formula <- model$formula$formula
  outcome_var <- all.vars(stats::formula(full_formula))[1]
  
  # Extract se() part
  se_expr <- full_formula[[2]]
  se_term <- deparse(se_expr[[3]][[2]])
  
  # Reconstruct formula with updated grouping
  formula_str <- glue::glue("{outcome_var} | se({se_term}) ~ 1 + (1 | {studyvar_name})")
  new_formula <- stats::as.formula(formula_str)
  
  stats::update(
    object = model,
    formula = new_formula,
    newdata = data,
    recompile = FALSE
  )
}

#' @export
update_model.meta_stan <- function(model, data, ...) {
  # MetaStan doesn't support updating
  # We return the original model; the adapter maps studies by position
  rlang::warn(
    c(
      "MetaStan models cannot be updated with new data.",
      "i" = "Study mapping uses row position from original data.",
      "i" = "Ensure data row order matches the order used when fitting the model."
    ),
    .frequency = "once",
    .frequency_id = "metastan_no_update"
  )
  model
}

# =============================================================================
# REFIT MODEL FOR SUBGROUPS (S3 Generic)
# =============================================================================

#' Refit Model for Subgroup Analysis
#'
#' @description
#' S3 generic for refitting a model on a subset of data (for subgroup analysis).
#' Each model type implements its own refitting logic.
#'
#' @param model A fitted model object
#' @param data Subset of data for this subgroup
#' @param original_model The original full model (for reference/priors)
#' @param ... Additional arguments passed to methods
#'
#' @return A refitted model object, or NULL if refitting not supported
#'
#' @keywords internal
#' @noRd
refit_model_subgroup <- function(model, data, original_model = NULL, ...) {
  
  UseMethod("refit_model_subgroup")
}

#' @export
refit_model_subgroup.default <- function(model, data, original_model = NULL, ...) {
  # Default: return NULL to indicate refitting not supported
  NULL
}

#' @export
refit_model_subgroup.brmsfit <- function(model, data, original_model = NULL, ...) {
  # brms supports direct model updating
  stats::update(model, newdata = data, recompile = FALSE, refresh = 0)
}

#' @export
refit_model_subgroup.meta_stan <- function(model, data, original_model = NULL, ...) {
  # MetaStan requires rebuilding the data object and refitting
  
  # Check if data has the required columns for create_MetaStan_dat
  # We need to detect whether this is arm-level or contrast-level data
  
  has_arm_cols <- all(c("r1", "r2", "n1", "n2") %in% names(data)) ||
    all(c("responders1", "responders2", "sampleSize1", "sampleSize2") %in% names(data))
  
  has_contrast_cols <- all(c("y", "se") %in% names(data)) ||
    all(c("yi", "sei") %in% names(data))
  
  if (!has_arm_cols && !has_contrast_cols) {
    rlang::warn(
      c(
        "Cannot refit MetaStan model for subgroup: data format not recognized.",
        "i" = "Data must contain either arm-level (r1, r2, n1, n2) or contrast-level (y, se) columns.",
        "i" = "Returning NULL - subgroup will use draws from overall model."
      ),
      .frequency = "once",
      .frequency_id = "metastan_subgroup_data"
    )
    return(NULL)
  }
  
  tryCatch({
    if (has_arm_cols) {
      # Detect column naming convention
      if (all(c("r1", "r2", "n1", "n2") %in% names(data))) {
        arm_vars <- c(responders = "r", sampleSize = "n")
      } else {
        arm_vars <- c(responders = "responders", sampleSize = "sampleSize")
      }
      
      # Create new MetaStan data object
      dat_subgroup <- MetaStan::create_MetaStan_dat(
        dat = as.data.frame(data),
        armVars = arm_vars
      )
      
      # Refit with same settings as original
      # Extract settings from original model if possible
      fit_subgroup <- MetaStan::meta_stan(
        data = dat_subgroup,
        likelihood = "binomial",
        mu_prior = c(0, 10),
        theta_prior = c(0, 100),
        tau_prior = 0.5,
        tau_prior_dist = "half-normal",
        chains = 4,
        iter = 2000,
        warmup = 1000,
        adapt_delta = 0.95
      )
      
      return(fit_subgroup)
    }
    # Add contrast-level handling if needed in future
    
  }, error = function(e) {
    rlang::warn(
      c(
        "Failed to refit MetaStan model for subgroup.",
        "i" = paste0("Error: ", e$message),
        "i" = "Returning NULL - subgroup will use draws from overall model."
      )
    )
    return(NULL)
  })
}

#' @export
refit_model_subgroup.bayesmeta <- function(model, data, original_model = NULL, ...) {
  # bayesmeta can be refit using the bayesmeta() function
  
  # Check for required columns
  if (!all(c("y", "se") %in% names(data)) && !all(c("yi", "sei") %in% names(data))) {
    rlang::warn(
      c(
        "Cannot refit bayesmeta model for subgroup: missing y/se columns.",
        "i" = "Returning NULL - subgroup will use draws from overall model."
      )
    )
    return(NULL)
  }
  
  tryCatch({
    # Get effect and SE columns
    y_col <- if ("y" %in% names(data)) "y" else "yi"
    se_col <- if ("se" %in% names(data)) "se" else "sei"
    
    # Refit bayesmeta model
    fit_subgroup <- bayesmeta::bayesmeta(
      y = data[[y_col]],
      sigma = data[[se_col]],
      labels = data$Author
    )
    
    return(fit_subgroup)
    
  }, error = function(e) {
    rlang::warn(
      c(
        "Failed to refit bayesmeta model for subgroup.",
        "i" = paste0("Error: ", e$message)
      )
    )
    return(NULL)
  })
}

#' @export
refit_model_subgroup.stanreg <- function(model, data, original_model = NULL, ...) {
  # rstanarm supports updating with new data
  
  tryCatch({
    fit_subgroup <- stats::update(model, data = data, refresh = 0)
    return(fit_subgroup)
    
  }, error = function(e) {
    rlang::warn(
      c(
        "Failed to refit rstanarm model for subgroup.",
        "i" = paste0("Error: ", e$message),
        "i" = "Returning NULL - subgroup will use draws from overall model."
      )
    )
    return(NULL)
  })
}

#' @export
refit_model_subgroup.stanfit <- function(model, data, original_model = NULL, ...) {
  # Raw rstan models cannot be easily refit without access to the original
  
  # Stan code and data preparation. Return NULL to use fallback.
  
  rlang::warn(
    c(
      "Cannot automatically refit rstan (stanfit) models for subgroup analysis.",
      "i" = "Subgroup pooled effects will be approximated from overall model draws.",
      "i" = "For proper subgroup analysis, consider fitting separate models per subgroup,",
      "i" = "or use brms/rstanarm which support automatic refitting."
    ),
    .frequency = "once",
    .frequency_id = "rstan_subgroup_refit"
  )
  
  return(NULL)
}

#' @export
refit_model_subgroup.CmdStanMCMC <- function(model, data, original_model = NULL, ...) {
  # CmdStanR models cannot be easily refit without access to the original
  # Stan code and data preparation. Return NULL to use fallback.
  
  rlang::warn(
    c(
      "Cannot automatically refit cmdstanr (CmdStanMCMC) models for subgroup analysis.",
      "i" = "Subgroup pooled effects will be approximated from overall model draws.",
      "i" = "For proper subgroup analysis, consider fitting separate models per subgroup,",
      "i" = "or use brms/rstanarm which support automatic refitting."
    ),
    .frequency = "once",
    .frequency_id = "cmdstanr_subgroup_refit"
  )
  
  return(NULL)
}