#' Internal function to extract draws from the posterior (multi-model version)
#'
#' @description
#' This is the updated version of forest.data_fn that works with any
#' supported model type via the extract_draws() generic.
#'
#' @param data Data frame with study information (must have Author column)
#' @param model Fitted model object
#' @param measure Effect measure type (needed for yi/vi calculation)
#' @param subgroup Logical for subgroup analysis
#' @param sort_studies_by How to sort studies
#' @param subgroup_order Order of subgroups
#' @param add_pred Whether to add prediction draws
#' @param add_pred_subgroup Whether to add prediction draws per subgroup
#'
#' @noRd
forest.data_fn <- function(data,
                           model,
                           measure = "OR",
                           subgroup = FALSE,
                           sort_studies_by = "author",
                           subgroup_order = NULL,
                           add_pred = FALSE,
                           add_pred_subgroup = FALSE) {
  
  model_type <- detect_model_type(model)
  
  # Ensure yi and vi exist in data (calculate from raw data if needed)
  data <- ensure_yi_vi(data, measure = measure)
  
  if (subgroup == FALSE && "Subgroup" %in% names(data)) {
    data <- data |> dplyr::select(-Subgroup)
  }
  
  if (subgroup == FALSE) {
    # =========================================================================
    # NON-SUBGROUP WORKFLOW
    # =========================================================================
    
    # Use the unified extract_draws function
    effect.draws <- extract_draws(model, data, add_pred = add_pred)
    
    # Join with study metadata - use any_of for optional columns
    effect.draws <- effect.draws |>
      dplyr::ungroup() |>
      dplyr::left_join(
        dplyr::select(data, Author, dplyr::any_of(c("Author_original", "Year", "yi", "vi"))),
        by = "Author"
      ) |>
      sort_studies_fn(sort_studies_by)
    
    # Ensure yi and vi columns exist after join
    if (!"yi" %in% names(effect.draws)) effect.draws$yi <- NA_real_
    if (!"vi" %in% names(effect.draws)) effect.draws$vi <- NA_real_
    
  } else {
    # =========================================================================
    # SUBGROUP WORKFLOW
    # =========================================================================
    
    # Nest data by subgroup
    subgroup_df <- data |>
      dplyr::group_by(Subgroup) |>
      tidyr::nest() |>
      dplyr::mutate(
        study_count = purrr::map_int(data, nrow)
      )
    
    # Attempt to refit model for each subgroup
    subgroup_df <- subgroup_df |>
      dplyr::mutate(
        subgroup_model = purrr::map(data, function(subgroup_data) {
          if (nrow(subgroup_data) < 2) {
            # Cannot fit meta-analysis with < 2 studies
            return(NULL)
          }
          refit_model_subgroup(model, subgroup_data, original_model = model)
        }),
        refit_success = purrr::map_lgl(subgroup_model, ~ !is.null(.x))
      )
    
    # Check if any refits failed
    if (any(!subgroup_df$refit_success & subgroup_df$study_count >= 2)) {
      failed_subgroups <- subgroup_df |>
        dplyr::filter(!refit_success & study_count >= 2) |>
        dplyr::pull(Subgroup)
      
      rlang::inform(
        c(
          paste0("Could not refit model for subgroup(s): ", paste(failed_subgroups, collapse = ", ")),
          "i" = "These subgroups will use study-level draws from the overall model.",
          "i" = "Subgroup pooled effects will be approximated."
        )
      )
    }
    
    # Extract draws for each subgroup
    study.effect.draws <- subgroup_df |>
      dplyr::mutate(
        effect_draws = purrr::pmap(
          list(subgroup_model, data, study_count, refit_success),
          function(sub_model, sub_data, n_studies, success) {
            
            if (success && !is.null(sub_model)) {
              # Model was successfully refit - extract draws from subgroup model
              combined <- extract_draws(sub_model, sub_data, 
                                        add_pred = (add_pred && add_pred_subgroup && n_studies > 1))
            } else {
              # Fallback: extract study draws from overall model, compute pooled manually
              all_draws <- extract_draws(model, data, add_pred = FALSE)
              
              # Filter to studies in this subgroup
              combined <- all_draws |>
                dplyr::filter(Author %in% sub_data$Author)
              
              if (n_studies >= 2) {
                # Compute approximate pooled effect for this subgroup
                subgroup_pooled <- combined |>
                  dplyr::group_by(.chain, .iteration, .draw) |>
                  dplyr::summarise(
                    b_Intercept = mean(b_Intercept),
                    sd_Author__Intercept = mean(sd_Author__Intercept, na.rm = TRUE),
                    .groups = "drop"
                  ) |>
                  dplyr::mutate(
                    Author = "Pooled Effect",
                    r_Author = NA_real_
                  )
                
                combined <- dplyr::bind_rows(combined, subgroup_pooled)
              }
            }
            
            # Handle single-study subgroups
            if (n_studies == 1) {
              combined <- combined |>
                dplyr::mutate(Author = dplyr::if_else(
                  Author == "Pooled Effect", 
                  "No Pooled Effect", 
                  Author
                ))
            }
            
            # Join metadata and sort - use any_of for optional columns
            combined |>
              dplyr::ungroup() |>
              dplyr::left_join(
                dplyr::select(sub_data, Author, dplyr::any_of(c("Author_original", "Year", "yi", "vi"))),
                by = "Author"
              ) |>
              sort_studies_fn(sort_studies_by)
          }
        )
      ) |>
      tidyr::unnest(effect_draws) |>
      dplyr::select(-data, -subgroup_model, -study_count, -refit_success)
    
    # Overall effect draws
    overall.effect.draws <- extract_draws(model, data, add_pred = FALSE) |>
      dplyr::filter(Author == "Pooled Effect") |>
      dplyr::mutate(
        Author = "Overall Effect",
        Author_original = NA_character_,
        Subgroup = "Overall",
        Year = NA_character_,
        yi = NA_real_,
        vi = NA_real_
      )
    
    effect.draws <- dplyr::bind_rows(study.effect.draws, overall.effect.draws)
    
    # Add overall prediction if requested
    if (isTRUE(add_pred)) {
      overall.pred.draws <- extract_draws(model, data, add_pred = TRUE) |>
        dplyr::filter(Author == "Prediction") |>
        dplyr::mutate(
          Author_original = NA_character_,
          Subgroup = "Overall",
          Year = NA_character_,
          yi = NA_real_,
          vi = NA_real_
        )
      
      effect.draws <- dplyr::bind_rows(effect.draws, overall.pred.draws)
    }
    
    # Apply subgroup ordering
    if (!is.null(subgroup_order)) {
      effect.draws <- effect.draws |>
        dplyr::mutate(Subgroup = factor(Subgroup, levels = subgroup_order)) |>
        dplyr::arrange(Subgroup) |>
        dplyr::mutate(
          Subgroup = dplyr::case_when(
            is.na(Subgroup) & Author == "Overall Effect" ~ "Overall",
            .default = as.character(Subgroup)
          )
        )
    }
    
    # Ensure yi and vi columns exist after all operations
    if (!"yi" %in% names(effect.draws)) effect.draws$yi <- NA_real_
    if (!"vi" %in% names(effect.draws)) effect.draws$vi <- NA_real_
  }
  
  return(effect.draws)
}

#' Internal function to summarise data for forest plot
#'
#' @noRd
forest.data.summary_fn <- function(spread_df,
                                   data,
                                   measure,
                                   sort_studies_by = "author",
                                   subgroup = FALSE,
                                   add_pred = FALSE,
                                   add_pred_subgroup = FALSE) {
  # Get effect size properties
  props <- get_measure_properties(measure)
  
  if (isFALSE(subgroup)){
    # Study Summaries
    forest.data <- spread_df |>
      dplyr::group_by(Author) |>
      tidybayes::median_qi(b_Intercept)
    
    # Tau Summary
    tau.summary <- spread_df |>
      dplyr::group_by(Author) |>
      tidybayes::median_qi(sd_Author__Intercept)
    
  } else if (subgroup == TRUE) {
    forest.data <- spread_df |>
      dplyr::group_by(Subgroup, Author) |>
      tidybayes::median_qi(b_Intercept)
    
    # Tau Summary
    tau.summary <- spread_df |>
      dplyr::group_by(Subgroup, Author) |>
      tidybayes::median_qi(sd_Author__Intercept)
  }
  
  # Select only needed join columns — include Author_original so it survives
  # Use any_of() for all optional columns (yi, vi, RoB columns)
  join_vars <- unique(c("Author", "Author_original", "Year", props$data_cols))
  rob_cols <- c("D1", "D2", "D3", "D4", "D5", "Overall")
  optional_cols <- c("yi", "vi", rob_cols)
  
  forest.data.summary <- forest.data |>
    dplyr::left_join(
      data |> dplyr::select(dplyr::any_of(c(join_vars, optional_cols))), 
      by = "Author"
    ) |>
    dplyr::left_join(tau.summary, by = "Author", suffix = c("", "_sd")) |> 
    sort_studies_fn(sort_studies_by)
  
  # Ensure yi and vi columns exist (set to NA if missing)
  if (!"yi" %in% names(forest.data.summary)) {
    forest.data.summary$yi <- NA_real_
  }
  if (!"vi" %in% names(forest.data.summary)) {
    forest.data.summary$vi <- NA_real_
  }
  
  # Ensure RoB columns exist (set to NA if missing)
  for (col in rob_cols) {
    if (!col %in% names(forest.data.summary)) {
      forest.data.summary[[col]] <- NA_character_
    }
  }
  
  # Add formatted effect estimates
  forest.data.summary <- forest.data.summary |>
    dplyr::mutate(
      weighted_effect = if (measure %in% c("MD", "SMD")) {
        paste0(sprintf("%.2f", b_Intercept), " [", sprintf("%.2f", .lower), ", ", sprintf("%.2f", .upper), "]")
      } else {
        paste0(sprintf("%.2f", exp(b_Intercept)), " [", sprintf("%.2f", exp(.lower)), ", ", sprintf("%.2f", exp(.upper)), "]")
      },
      # Handle NA yi/vi gracefully
      unweighted_effect = dplyr::case_when(
        is.na(yi) | is.na(vi) ~ NA_character_,
        measure %in% c("MD", "SMD") ~ paste0(
          sprintf("%.2f", yi), " [", 
          sprintf("%.2f", yi - 1.96 * sqrt(vi)), ", ", 
          sprintf("%.2f", yi + 1.96 * sqrt(vi)), "]"
        ),
        TRUE ~ paste0(
          sprintf("%.2f", exp(yi)), " [", 
          sprintf("%.2f", exp(yi - 1.96 * sqrt(vi))), ", ", 
          sprintf("%.2f", exp(yi + 1.96 * sqrt(vi))), "]"
        )
      ),
      # Replace NA with tau summary for pooled rows
      unweighted_effect = dplyr::case_when(
        is.na(unweighted_effect) & Author %in% c("Pooled Effect", "Overall Effect") ~
          paste0("\u03c4 = ", sprintf("%.2f", sd_Author__Intercept),
                 " [", sprintf("%.2f", .lower_sd), ", ", sprintf("%.2f", .upper_sd), "]"),
        is.na(unweighted_effect) ~ "",
        TRUE ~ unweighted_effect
      )
    )
  
  # For the Prediction row, blank out the unweighted_effect (no tau to display)
  forest.data.summary <- forest.data.summary |>
    dplyr::mutate(
      unweighted_effect = dplyr::if_else(
        as.character(Author) == "Prediction",
        "",
        unweighted_effect
      )
    )
  
  # Add study group summary columns depending on measure type
  if (measure %in% c("MD", "SMD")) {
    # Check if required columns exist
    has_cont_cols <- all(c("N_Intervention", "N_Control", "Mean_Intervention", 
                           "Mean_Control", "SD_Intervention", "SD_Control") %in% 
                           names(forest.data.summary))
    
    if (has_cont_cols) {
      forest.data.summary <- forest.data.summary |>
        dplyr::mutate(
          N_int = dplyr::case_when(
            Author == "Pooled Effect" ~ as.character(sum(N_Intervention, na.rm = TRUE)),
            TRUE ~ as.character(N_Intervention)),
          int_mean_sd = dplyr::case_when(
            Author %in% c("Pooled Effect", "Prediction") ~ NA_character_,
            TRUE ~ paste0(sprintf("%.2f", Mean_Intervention), " (", sprintf("%.2f", SD_Intervention), ")")),
          N_ctrl = dplyr::case_when(
            Author == "Pooled Effect" ~ as.character(sum(N_Control, na.rm = TRUE)),
            TRUE ~ as.character(N_Control)),
          ctrl_mean_sd = dplyr::case_when(
            Author %in% c("Pooled Effect", "Prediction") ~ NA_character_,
            TRUE ~ paste0(sprintf("%.2f", Mean_Control), " (", sprintf("%.2f", SD_Control), ")"))
        )
    } else {
      forest.data.summary <- forest.data.summary |>
        dplyr::mutate(
          N_int = NA_character_,
          int_mean_sd = NA_character_,
          N_ctrl = NA_character_,
          ctrl_mean_sd = NA_character_
        )
    }
  } else {
    # Check if event/n columns exist before using them
    has_event_cols <- all(c("Event_Control", "N_Control", "Event_Intervention", "N_Intervention") %in% 
                            names(forest.data.summary))
    
    if (has_event_cols) {
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
    } else {
      # Create empty columns if data doesn't have events
      forest.data.summary <- forest.data.summary |>
        dplyr::mutate(
          control_outcome_frac = NA_character_,
          int_outcome_frac = NA_character_
        )
    }
  }
  
  # Remove rows where Author == "No Pooled Effect"
  forest.data.summary <- forest.data.summary |>
    dplyr::filter(Author != "No Pooled Effect")
  
  return(forest.data.summary)
}