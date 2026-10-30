#' Internal function to extract draws from the posterior
#'
#' @description
#' Uses `extract_forest_draws()` from the adapter layer to extract posterior
#' draws from any supported model backend. Returns draws in canonical format
#' with columns: .draw, Author, theta_study, eta, tau.
#'
#' @param data Data frame with study information (must contain Author column).
#' @param model A fitted model object (any supported backend).
#' @param subgroup Logical. Whether to perform subgroup analysis.
#' @param sort_studies_by How to sort studies.
#' @param subgroup_order Character vector of subgroup ordering.
#' @param effect_type Character: "random_effect" or "common_effect".
#' @param re_min_k Passed to subgroup refits (bayesma models only).
#'
#' @noRd
forest.data_fn <- function(data,
                           model,
                           subgroup = FALSE,
                           sort_studies_by = "author",
                           subgroup_order = NULL,
                           effect_type = "random_effect",
                           add_pred = FALSE,
                           add_pred_subgroup = FALSE,
                           re_min_k = NULL) {

  if (subgroup == FALSE && "Subgroup" %in% names(data)) {
    data <- data |> dplyr::select(-Subgroup)
  }

  if (subgroup == FALSE) {
    # No subgroup: extract draws using adapter
    draws <- extract_forest_draws(model, data)
    if (isTRUE(add_pred)) draws <- add_prediction_draws(draws)

    effect.draws <- draws |>
      dplyr::mutate(Author = stringr::str_replace_all(Author, "\\.", " ")) |>
      dplyr::ungroup() |>
      dplyr::left_join(
        dplyr::select(data, Author, dplyr::any_of(c("Author_original", "Year", "yi", "vi"))),
        by = "Author") |>
      sort_studies_fn(sort_studies_by)

  } else {
    # With subgroup
    full_data <- data  # capture before nesting shadows the name
    subgroup_df <- data |>
      dplyr::group_by(Subgroup) |>
      tidyr::nest() |>
      dplyr::mutate(
        subgroup_model = purrr::map(data, ~ refit_for_subgroup(model, .x, original_data = full_data, re_min_k = re_min_k)),
        study_count = purrr::map_int(data, nrow)
      )

    # Overall effect draws
    overall_draws <- extract_forest_draws(model, data)
    if (isTRUE(add_pred)) overall_draws <- add_prediction_draws(overall_draws)

    overall.effect.draws <- overall_draws |>
      dplyr::filter(Author %in% c("Pooled Effect", "Prediction")) |>
      dplyr::mutate(
        Author    = dplyr::if_else(Author == "Pooled Effect", "Overall Effect", Author),
        Subgroup  = "Overall",
        Author_original = NA_character_,
        Year      = NA_integer_,
        yi        = NA_real_,
        vi        = NA_real_
      )

    # Extract draws for each subgroup model
    study.effect.draws <- subgroup_df |>
      dplyr::mutate(
        effect_draws = purrr::pmap(
          list(subgroup_model, data, study_count, Subgroup),   # <- add Subgroup
          function(sg_model, sg_data, n_studies, sg_label) {   # <- receive it as sg_label
            # If refit returned NULL, filter from overall model draws
            if (is.null(sg_model)) {
              sg_study_names <- as.character(sg_data$Author)

              study_draws <- overall_draws |>
                dplyr::filter(Author %in% sg_study_names)

              # Approximate subgroup pooled from study draws
              pooled_label <- ifelse(n_studies == 1, "No Pooled Effect", "Pooled Effect")
              pooled_draws <- study_draws |>
                dplyr::group_by(.draw) |>
                dplyr::summarise(
                  theta_study = mean(theta_study),
                  eta         = NA_real_,
                  tau         = mean(tau, na.rm = TRUE),
                  .groups     = "drop"
                ) |>
                dplyr::mutate(Author = pooled_label)

              combined <- dplyr::bind_rows(study_draws, pooled_draws)
            } else {
              combined <- extract_forest_draws(sg_model, sg_data)

              # Rename pooled label if only 1 study
              pooled_label <- ifelse(n_studies == 1, "No Pooled Effect", "Pooled Effect")
              if (pooled_label != "Pooled Effect") {
                combined <- combined |>
                  dplyr::mutate(Author = dplyr::if_else(
                    Author == "Pooled Effect", pooled_label, Author
                  ))
              }
            }

            if (isTRUE(add_pred_subgroup)) combined <- add_prediction_draws(combined)

            combined |>
              dplyr::mutate(Author = stringr::str_replace_all(Author, "\\.", " ")) |>
              dplyr::ungroup() |>
              dplyr::left_join(
                dplyr::select(sg_data, Author, dplyr::any_of(c("Author_original", "Year", "yi", "vi"))),
                by = "Author"
              ) |>
              sort_studies_fn(sort_studies_by)
          })
      ) |>
      tidyr::unnest(effect_draws) |>
      dplyr::select(-data, -subgroup_model, -study_count)

    effect.draws <- dplyr::bind_rows(study.effect.draws, overall.effect.draws)

    # Custom group order for Subgroup column
    if (!is.null(subgroup_order)) {
      effect.draws <- effect.draws |>
        dplyr::mutate(Subgroup = factor(Subgroup, levels = subgroup_order)) |>
        dplyr::arrange(Subgroup) |>
        dplyr::mutate(Subgroup = dplyr::case_when(
          is.na(Subgroup) & Author == "Overall Effect" ~ "Overall",
          .default = Subgroup
        ))
    }
  }
  return(effect.draws)
}

#' Internal function to summarise data for forest plot
#'
#' @param spread_df Draws tibble with canonical columns.
#' @param data Original data frame.
#' @param measure Effect measure.
#' @param sort_studies_by Sorting method.
#' @param subgroup Logical.
#' @param effect_type "random_effect" or "common_effect".
#' @param rob_cols Risk of bias columns to carry through to the summary.
#'
#' @noRd
forest.data.summary_fn <- function(spread_df,
                                   data,
                                   measure,
                                   sort_studies_by = "author",
                                   subgroup = FALSE,
                                   effect_type = "random_effect",
                                   rob_cols = character()) {
  # Get effect size properties
  props <- get_measure_properties(measure)

  if (isFALSE(subgroup)) {
    # Study Summaries
    forest.data <- spread_df |>
      dplyr::group_by(Author) |>
      tidybayes::median_qi(theta_study)

    # Tau Summary (only for RE models)
    if (effect_type == "random_effect") {
      tau.summary <- spread_df |>
        dplyr::group_by(Author) |>
        tidybayes::median_qi(tau)
    }

  } else if (subgroup == TRUE) {
    forest.data <- spread_df |>
      dplyr::group_by(Subgroup, Author) |>
      tidybayes::median_qi(theta_study)

    if (effect_type == "random_effect") {
      tau.summary <- spread_df |>
        dplyr::group_by(Subgroup, Author) |>
        tidybayes::median_qi(tau)
    }
  }

  # Select only needed join columns
  join_vars <- unique(c("Author", "Author_original", "Year", props$data_cols))
  forest.data.summary <- forest.data |>
    dplyr::left_join(data |> dplyr::select(dplyr::any_of(join_vars), yi, vi, dplyr::any_of(rob_cols)), by = "Author")

  if (effect_type == "random_effect") {
    forest.data.summary <- forest.data.summary |>
      dplyr::left_join(tau.summary, by = "Author", suffix = c("", "_tau"))
  }

  forest.data.summary <- forest.data.summary |>
    sort_studies_fn(sort_studies_by)

  # Add formatted effect estimates
  forest.data.summary <- forest.data.summary |>
    dplyr::mutate(
      weighted_effect = if (measure %in% c("MD", "SMD")) {
        paste0(sprintf("%.2f", theta_study), " [", sprintf("%.2f", .lower), ", ", sprintf("%.2f", .upper), "]")
      } else {
        paste0(sprintf("%.2f", exp(theta_study)), " [", sprintf("%.2f", exp(.lower)), ", ", sprintf("%.2f", exp(.upper)), "]")
      })

  # Unweighted effect (observed) and tau display
  if (effect_type == "random_effect") {
    forest.data.summary <- forest.data.summary |>
      dplyr::mutate(
        unweighted_effect = if (measure %in% c("MD", "SMD")) {
          paste0(sprintf("%.2f", yi), " [", sprintf("%.2f", yi - 1.96 * sqrt(vi)), ", ", sprintf("%.2f", yi + 1.96 * sqrt(vi)), "]")
        } else {
          paste0(sprintf("%.2f", exp(yi)), " [", sprintf("%.2f", exp(yi - 1.96 * sqrt(vi))), ", ", sprintf("%.2f", exp(yi + 1.96 * sqrt(vi))), "]")
        },
        unweighted_effect = ifelse(
          unweighted_effect == "NA [NA, NA]",
          paste0("\u03c4 = ", sprintf("%.2f", tau),
                 " [", sprintf("%.2f", .lower_tau), ", ", sprintf("%.2f", .upper_tau), "]"),
          unweighted_effect
        ),
        # Blank out unweighted_effect for Prediction rows (no tau to display)
        unweighted_effect = dplyr::if_else(
          as.character(Author) == "Prediction", "", unweighted_effect
        ))
  } else {
    # Common-effects: no tau, just observed effects (or blank for pooled)
    forest.data.summary <- forest.data.summary |>
      dplyr::mutate(
        unweighted_effect = if (measure %in% c("MD", "SMD")) {
          paste0(sprintf("%.2f", yi), " [", sprintf("%.2f", yi - 1.96 * sqrt(vi)), ", ", sprintf("%.2f", yi + 1.96 * sqrt(vi)), "]")
        } else {
          paste0(sprintf("%.2f", exp(yi)), " [", sprintf("%.2f", exp(yi - 1.96 * sqrt(vi))), ", ", sprintf("%.2f", exp(yi + 1.96 * sqrt(vi))), "]")
        },
        unweighted_effect = ifelse(
          unweighted_effect == "NA [NA, NA]", "\u2014", unweighted_effect
        ),
        unweighted_effect = dplyr::if_else(
          as.character(Author) == "Prediction", "", unweighted_effect
        ))
  }

  # Add study group summary columns depending on measure type
  if (measure %in% c("MD", "SMD")) {
    forest.data.summary <- forest.data.summary |>
      dplyr::mutate(
        N_int = dplyr::case_when(
          Author == "Pooled Effect" ~ as.character(sum(N_Intervention, na.rm = TRUE)),
          Author == "Prediction" ~ "",
          TRUE ~ as.character(N_Intervention)),
        int_mean_sd = dplyr::case_when(
          Author %in% c("Pooled Effect", "Prediction") ~ NA_character_,
          TRUE ~ paste0(sprintf("%.2f", Mean_Intervention), " (", sprintf("%.2f", SD_Intervention), ")")),
        N_ctrl = dplyr::case_when(
          Author == "Pooled Effect" ~ as.character(sum(N_Control, na.rm = TRUE)),
          Author == "Prediction" ~ "",
          TRUE ~ as.character(N_Control)),
        ctrl_mean_sd = dplyr::case_when(
          Author %in% c("Pooled Effect", "Prediction") ~ NA_character_,
          TRUE ~ paste0(sprintf("%.2f", Mean_Control), " (", sprintf("%.2f", SD_Control), ")"))
      )
  } else {
    forest.data.summary <- forest.data.summary |>
      dplyr::mutate(
        control_outcome_frac = dplyr::case_when(
          Author == "Pooled Effect" ~ paste0(
            sum(Event_Control[Author != "Pooled Effect"], na.rm = TRUE), "/",
            sum(N_Control[Author != "Pooled Effect"], na.rm = TRUE)),
          Author == "Prediction" ~ "",
          TRUE ~ paste0(Event_Control, "/", N_Control)),
        int_outcome_frac = dplyr::case_when(
          Author == "Pooled Effect" ~ paste0(
            sum(Event_Intervention[Author != "Pooled Effect"], na.rm = TRUE), "/",
            sum(N_Intervention[Author != "Pooled Effect"], na.rm = TRUE)),
          Author == "Prediction" ~ "",
          TRUE ~ paste0(Event_Intervention, "/", N_Intervention)))
  }

  # Remove rows where Author == "No Pooled Effect"
  forest.data.summary <- forest.data.summary |>
    dplyr::filter(Author != "No Pooled Effect")

  return(forest.data.summary)
}
