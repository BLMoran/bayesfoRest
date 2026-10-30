
#' Internal function to sort studies
#'
#' @noRd
sort_studies_fn <- function(.data, sort_studies_by = NULL) {
  sort_by <- if (is.null(sort_studies_by)) "author" else sort_studies_by
  has_subgroup <- "Subgroup" %in% names(.data) &&
    !all(is.na(.data$Subgroup)) &&
    dplyr::n_distinct(.data$Subgroup) > 1


  # Separate pooled effect and prediction rows
  special_rows <- c("Pooled Effect", "Prediction")
  pooled     <- .data |> dplyr::filter(Author == "Pooled Effect")
  prediction <- .data |> dplyr::filter(Author == "Prediction")
  studies    <- .data |> dplyr::filter(!Author %in% special_rows)

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
#'    This is needed so that ggplot2 maps each study to its own row on the
#'    y-axis.
#'
#' The original Author name is preserved in a new column `Author_original` so
#' it can be restored for display purposes later.
#'
#' @param data A data frame containing at least `Author` and `Year` columns.
#'
#' @return The data frame with unique `Author` values, disambiguated `Year`
#'   values where needed, and a new `Author_original` column containing the
#'   original (possibly duplicated) names.
#'
#' @noRd
make_authors_unique <- function(data, studyvar = "Author", year = NULL) {
  # resolve column names (accept bare or string)
  study_col <- rlang::as_name(rlang::ensym(studyvar))
  year_col  <- if (!is.null(year)) rlang::as_name(rlang::ensym(year)) else NULL

  if (!study_col %in% names(data)) {
    cli::cli_abort("Column {.field {study_col}} not found in data.")
  }

  # preserve original label for display
  data[["Author_original"]] <- data[[study_col]]

  # --- Step 1: Year disambiguation (only if a valid year column is present) ---
  if (!is.null(year_col) && year_col %in% names(data)) {
    data[[year_col]] <- as.character(data[[year_col]])
    combos <- paste(data[[study_col]], data[[year_col]], sep = "|||")
    dup <- combos[duplicated(combos)]
    for (combo in unique(dup)) {
      idx <- which(combos == combo)
      yr  <- strsplit(combo, "|||", fixed = TRUE)[[1]][2]
      data[[year_col]][idx] <- paste0(yr, letters[seq_along(idx)])
    }
  }

  # --- Step 2: Study-label disambiguation (always) ---
  dup_ids <- data[[study_col]][duplicated(data[[study_col]])]
  for (nm in unique(dup_ids)) {
    idx <- which(data[[study_col]] == nm)
    data[[study_col]][idx] <- paste0(nm, "_", letters[seq_along(idx)])
  }

  data
}


#' Enforce a minimum number of studies for random-effects models
#'
#' If `re_min_k` is set and the number of unique studies is below it,
#' a random-effects model_type is downgraded to common_effect with a warning.
#' Other model types are passed through unchanged.
#'
#' @noRd
enforce_re_min_k <- function(model_type, re_min_k, data, studyvar) {
  if (is.null(re_min_k)) return(model_type)

  if (!is.numeric(re_min_k) || length(re_min_k) != 1 || re_min_k < 1) {
    cli::cli_abort(
      "{.arg re_min_k} must be a single positive number.",
      call = rlang::caller_env()
    )
  }

  resolved <- if (length(model_type) > 1) model_type[1] else model_type
  if (resolved != "random_effect") return(model_type)

  k <- length(unique(data[[studyvar]]))
  if (k >= re_min_k) return(model_type)

  cli::cli_warn(c(
    "!" = "Switching to {.val common_effect}: only {k} stud{?y/ies} \\
           (< {.arg re_min_k} = {re_min_k}).",
    "i" = "Random-effects estimation is unstable below {re_min_k} studies."
  ))
  "common_effect"
}
