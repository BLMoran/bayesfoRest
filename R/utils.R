
#' Internal function to sort studies
#'
#' @noRd
sort_studies_fn <- function(.data, sort_studies_by = NULL) {
  sort_by <- if (is.null(sort_studies_by)) "author" else sort_studies_by
  has_subgroup <- "Subgroup" %in% names(.data) &&
    !all(is.na(.data$Subgroup)) &&
    dplyr::n_distinct(.data$Subgroup) > 1
  
  
  # Separate pooled effect if present
  pooled <- .data |> dplyr::filter(Author == "Pooled Effect")
  studies <- .data |> dplyr::filter(Author != "Pooled Effect")
  
  # Arrange studies
  studies <- switch(
    sort_by,
    "author" = studies |> dplyr::arrange(dplyr::across(any_of("Subgroup")), Author),
    "year" = studies |> dplyr::arrange(dplyr::across(any_of("Subgroup")), Year, Author),
    "effect" = studies |> dplyr::arrange(dplyr::across(any_of("Subgroup")), yi, Author),
    stop("Invalid value for `sort_studies_by`. Must be one of 'author', 'year', or 'effect'.")
  )
  
  # Recombine and set factor levels
  full <- dplyr::bind_rows(studies, pooled)
  
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


#' Internal function to update model
#'
#' @noRd
update_model <- function(model, data, studyvar) {
  studyvar_name <- rlang::as_name(rlang::ensym(studyvar))
  
  # Extract response variable from formula that includes se()
  full_formula <- model$formula$formula
  outcome_var <- all.vars(stats::formula(full_formula))[1]
  
  # Extract se() part
  se_expr <- full_formula[[2]]
  se_term <- deparse(se_expr[[3]][[2]])
  
  # Reconstruct new formula string with updated grouping variable
  formula_str <- glue::glue("{outcome_var} | se({se_term}) ~ 1 + (1 | {studyvar_name})")
  
  # Convert to formula
  new_formula <- stats::as.formula(formula_str)
  
  # Update model
  model <- stats::update(
    object = model,
    formula = new_formula,
    newdata = data,
    recompile = FALSE
  )
  return(model)
}