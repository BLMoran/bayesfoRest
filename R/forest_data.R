#' Internal function to extract draws from the posterior
#'
#' @noRd
forest.data_fn <- function(data,
                           model,
                           subgroup = FALSE,
                           sort_studies_by = "author",
                           subgroup_order = NULL,
                           add_pred = FALSE,
                           add_pred_subgroup = FALSE) {
  
  if (subgroup == FALSE && "Subgroup" %in% names(data)) {
    data <- data |> dplyr::select(-Subgroup)
  }
  if (subgroup == FALSE) {
    # No subgroup, evaluate as a whole
    study.draws <- tidybayes::spread_draws(model, r_Author[Author, ], b_Intercept) |>
      dplyr::mutate(b_Intercept = r_Author + b_Intercept)
    pooled.draws <- tidybayes::spread_draws(model, b_Intercept, sd_Author__Intercept) |>
      dplyr::mutate(Author = "Pooled Effect")
    
    effect.draws <- dplyr::bind_rows(study.draws, pooled.draws)
    
    # Generate prediction draws if requested
    if (isTRUE(add_pred)) {
      nd <- data.frame(Author = "new", sei = 0)
      pred_samples <- brms::posterior_predict(
        object = model,
        newdata = nd,
        re_formula = NULL,
        allow_new_levels = TRUE,
        sample_new_levels = "gaussian"
      )
      
      # Create prediction draws in the same structure as other draws
      pred.draws <- tidybayes::spread_draws(model, b_Intercept, sd_Author__Intercept) |>
        dplyr::mutate(
          Author = "Prediction",
          b_Intercept = as.vector(pred_samples)
        )
      effect.draws <- dplyr::bind_rows(effect.draws, pred.draws)
    }
    
    effect.draws <- effect.draws |>
      dplyr::mutate(Author = stringr::str_replace_all(Author, "\\.", " ")) |> 
      dplyr::ungroup() |>
      dplyr::left_join(dplyr::select(data, Author, Author_original, Year, yi, vi), by = "Author") |>
      sort_studies_fn(sort_studies_by)
  } else {
    # With subgroup
    subgroup_df <- data |>
      dplyr::group_by(Subgroup) |>
      tidyr::nest() |>
      dplyr::mutate(
        subgroup_model = purrr::map(data, ~ stats::update(model, newdata = .x)),
        # Calculate study count for each subgroup
        study_count = purrr::map_int(data, nrow)
      )
    
    # Create overall effect draws with matching structure to subgroup effect.draws
    overall.effect.draws <- tidybayes::spread_draws(model, b_Intercept, sd_Author__Intercept) |>
      dplyr::mutate(
        Author = "Overall Effect",
        Author_original = NA_character_,
        Subgroup = "Overall",
        r_Author = 0,           # No random effect for overall
        Year = NA_character_,   # Character to match disambiguated Year from make_authors_unique()
        yi = NA_real_,          # No individual effect size
        vi = NA_real_           # No individual variance
      )
    
    # Extract draws for each subgroup model
    study.effect.draws <- subgroup_df |>
      dplyr::mutate(
        effect_draws = purrr::pmap(list(subgroup_model, data, study_count), ~ {
          study <- tidybayes::spread_draws(..1, r_Author[Author, ], b_Intercept) |>
            dplyr::mutate(b_Intercept = r_Author + b_Intercept)
          
          # Use different label based on study count
          pooled_label <- ifelse(..3 == 1, "No Pooled Effect", "Pooled Effect")
          
          pooled <- tidybayes::spread_draws(..1, b_Intercept, sd_Author__Intercept) |>
            dplyr::mutate(Author = pooled_label)
          
          combined <- dplyr::bind_rows(study, pooled)
          
          # Generate prediction draws for this subgroup if requested
          if (isTRUE(add_pred) && isTRUE(add_pred_subgroup) && ..3 > 1) {
            nd <- data.frame(Author = "new", sei = 0)
            pred_samples <- brms::posterior_predict(
              object = ..1,
              newdata = nd,
              re_formula = NULL,
              allow_new_levels = TRUE,
              sample_new_levels = "gaussian"
            )
            
            pred <- tidybayes::spread_draws(..1, b_Intercept, sd_Author__Intercept) |>
              dplyr::mutate(
                Author = "Prediction",
                b_Intercept = as.vector(pred_samples)
              )
            combined <- dplyr::bind_rows(combined, pred)
          }
          
          combined <- combined |>
            dplyr::mutate(Author = stringr::str_replace_all(Author, "\\.", " ")) |> 
            dplyr::ungroup() |>
            dplyr::left_join(dplyr::select(..2, Author, Author_original, Year, yi, vi), by = "Author") |>
            sort_studies_fn(sort_studies_by)
          return(combined)
        })
      ) |>
      tidyr::unnest(effect_draws) |>
      dplyr::select(-data, -subgroup_model, -study_count)
    
    effect.draws <- dplyr::bind_rows(study.effect.draws, overall.effect.draws)
    
    # Add overall prediction draws when add_pred = TRUE
    if (isTRUE(add_pred)) {
      nd <- data.frame(Author = "new", sei = 0)
      pred_samples <- brms::posterior_predict(
        object = model,
        newdata = nd,
        re_formula = NULL,
        allow_new_levels = TRUE,
        sample_new_levels = "gaussian"
      )
      
      overall.pred.draws <- tidybayes::spread_draws(model, b_Intercept, sd_Author__Intercept) |>
        dplyr::mutate(
          Author = "Prediction",
          Author_original = NA_character_,
          Subgroup = "Overall",
          r_Author = 0,
          Year = NA_character_,
          yi = NA_real_,
          vi = NA_real_,
          b_Intercept = as.vector(pred_samples)
        )
      effect.draws <- dplyr::bind_rows(effect.draws, overall.pred.draws)
    }
    
    # Custom group order for Subgroup column
    if (!is.null(subgroup_order)) {
      effect.draws <- effect.draws |>
        dplyr::mutate(Subgroup = factor(Subgroup, levels = subgroup_order)) |>
        dplyr::arrange(Subgroup) |>
        dplyr::mutate(Subgroup =  dplyr::case_when(
          is.na(Subgroup) & Author == "Overall Effect" ~ "Overall",
          .default = Subgroup
        ))
    }
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
  join_vars <- unique(c("Author", "Author_original", "Year", props$data_cols))
  forest.data.summary <- forest.data |>
    dplyr::left_join(data |> dplyr::select(dplyr::any_of(join_vars), yi, vi, D1:Overall), by = "Author") |>
    dplyr::left_join(tau.summary, by = "Author", suffix = c("", "_sd")) |> 
    sort_studies_fn(sort_studies_by)
  
  # Add formatted effect estimates
  forest.data.summary <- forest.data.summary |>
    dplyr::mutate(
      weighted_effect = if (measure %in% c("MD", "SMD")) {
        paste0(sprintf("%.2f", b_Intercept), " [", sprintf("%.2f", .lower), ", ", sprintf("%.2f", .upper), "]")
      } else {
        paste0(sprintf("%.2f", exp(b_Intercept)), " [", sprintf("%.2f", exp(.lower)), ", ", sprintf("%.2f", exp(.upper)), "]")
      },
      unweighted_effect = if (measure %in% c("MD", "SMD")) {
        paste0(sprintf("%.2f", yi), " [", sprintf("%.2f", yi - 1.96 * sqrt(vi)), ", ", sprintf("%.2f", yi + 1.96 * sqrt(vi)), "]")
      } else {
        paste0(sprintf("%.2f", exp(yi)), " [", sprintf("%.2f", exp(yi - 1.96 * sqrt(vi))), ", ", sprintf("%.2f", exp(yi + 1.96 * sqrt(vi))), "]")
      },
      unweighted_effect = ifelse(unweighted_effect == "NA [NA, NA]",
                                 paste0("\u03c4 = ", sprintf("%.2f", forest.data.summary$sd_Author__Intercept),
                                        " [", sprintf("%.2f", forest.data.summary$.lower_sd), ", ", sprintf("%.2f", forest.data.summary$.upper_sd), "]"),
                                 unweighted_effect))
  
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

#' Extract Study-Level Effects
#'
#' @description
#' Extracts study-specific effect estimates combining fixed and random effects
#' from a Bayesian meta-analysis model.
#'
#' @param model A fitted brmsfit object
#' @param data Data frame containing study data
#' @param study_var Name of the study identifier variable
#' @param measure Effect measure type
#'
#' @return List of numeric vectors, one per study, containing posterior draws
#'
#' @keywords internal
#' @noRd
extract_study_level_effects <- function(model, data, study_var, measure) {
  
  # Get random effects
  ranef_draws <- brms::ranef(model, summary = FALSE)
  
  # Get the study-level effects
  if (study_var %in% names(ranef_draws)) {
    study_effects <- ranef_draws[[study_var]][, , "Intercept"]
  } else {
    # If no random effects, use fixed effect only
    study_effects <- matrix(0, nrow = nrow(as.matrix(model)), ncol = length(unique(data[[study_var]])))
    colnames(study_effects) <- unique(data[[study_var]])
  }
  
  # Get overall intercept
  intercept_draws <- tidybayes::spread_draws(model, b_Intercept)$b_Intercept
  
  # Combine for each study
  purrr::map(seq_len(nrow(data)), function(i) {
    study_name <- as.character(data[[study_var]][i])
    if (study_name %in% colnames(study_effects)) {
      draws <- intercept_draws + study_effects[, study_name]
    } else {
      draws <- intercept_draws
    }
    
    # Transform if necessary
    if (measure %in% c("OR", "RR", "HR", "IRR")) {
      exp(draws)
    } else {
      draws
    }
  })
}