#' Create Bayesian Forest Plot for Meta-Analysis
#'
#' @description
#' This function creates a Bayesian forest plot for meta-analysis
#' using data prepared with metafor and a fitted brms model. It supports various effect measures and can handle
#' both single and subgroup analyses with customizable visualization options.
#'
#' @param model A fitted brms model object (class 'brmsfit') containing the 
#'   Bayesian meta-analysis results.
#' @param data A data frame containing the study data used for the meta-analysis.
#' @param measure Character string specifying the effect measure. Must be one of:
#'   "OR" (Odds Ratio), "HR" (Hazard Ratio), "RR" (Risk Ratio), 
#'   "IRR" (Incidence Rate Ratio), "MD" (Mean Difference), or "SMD" (Standardized Mean Difference).
#' @param studyvar Column name containing study identifiers/authors. Default is NULL.
#' @param year Column name containing publication years. Default is NULL.
#' @param c_n Column name containing control group sample sizes. Required for OR, RR, MD, SMD.
#' @param i_n Column name containing intervention group sample sizes. Required for OR, RR, MD, SMD.
#' @param c_event Column name containing control group event counts. Required for OR, RR, IRR.
#' @param i_event Column name containing intervention group event counts. Required for OR, RR, IRR.
#' @param c_mean Column name containing control group means. Required for MD, SMD.
#' @param i_mean Column name containing intervention group means. Required for MD, SMD.
#' @param c_sd Column name containing control group standard deviations. Required for MD, SMD.
#' @param i_sd Column name containing intervention group standard deviations. Required for MD, SMD.
#' @param c_time Column name containing control group time periods. Required for IRR.
#' @param i_time Column name containing intervention group time periods. Required for IRR.
#' @param sort_studies_by Character string specifying how to sort studies. 
#'   Options: "author" (default), "year", or "effect".
#' @param subgroup Logical indicating whether to create subgroup analysis. Default is FALSE.
#' @param sort_subgroup_by Character string or vector specifying subgroup ordering.
#'   Options: "alphabetical" (default), "effect", or custom character vector of subgroup names.
#' @param label_outcome Character string for outcome label. Default is "Outcome".
#' @param label_control Character string for control group label. Default is "Control".
#' @param label_intervention Character string for intervention group label. Default is "Intervention".
#' @param title Character string for the plot title. Default is NULL (no title).
#' @param subtitle Character string for the plot subtitle. Default is NULL (no subtitle).
#' @param title_align Character string specifying title alignment. Options: "left" (default), "center"/"centre", "right".
#' @param xlim Numeric vector of length 2 specifying x-axis limits. Default is NULL (auto-scaled).
#' @param x_breaks Numeric vector specifying custom x-axis break points. Default is NULL (uses measure-specific defaults).
#' @param shrinkage_output Character string specifying shrinkage visualization. 
#'   Options: "density" (default) or "pointinterval".
#' @param null_value Numeric value specifying x-axis value of the null line. Default is NULL (uses measure-specific defaults).
#' @param add_rope Logical indicating whether to add ROPE (Region of Practical Equivalence). Default is FALSE.
#' @param rope_value Numeric vector of length 2 specifying ROPE range, or single value for symmetric range around null. 
#'   Default is NULL (uses Kruschke's recommendations: OR/HR/RR/IRR = c(0.9, 1.1), SMD = c(-0.1, 0.1), MD requires specification).
#' @param rope_color Color for ROPE shading. Default is transparent grey.
#' @param color_palette Character vector of colors for the plot. Default is NULL.
#' @param color_study_posterior_null_left Color for left side of null study posteriors. Default is "deepskyblue".
#' @param color_study_posterior_null_right Color for right side of null study posteriors. Default is "violet".
#' @param color_study_posterior Color for study posterior densities. Default is "dodgerblue".
#' @param color_study_posterior_outline Color for study posterior outlines. Default is "blue".
#' @param color_overall_posterior Color for overall posterior. Default is "blue".
#' @param color_pointinterval Color for point intervals. Default is "purple".
#' @param color_shrinkage_outline Color for shrinkage plot outlines. Default is "purple".
#' @param color_shrinkage_fill Color for shrinkage plot fill. Default is NULL.
#' @param plot_width Numeric value specifying the relative width of the plot component. Default is 4.
#' @param add_rob Logical indicating whether to add Risk of Bias assessment. Default is FALSE.
#' @param rob_tool Character string specifying RoB tool. Options: "rob2" (default).
#' @param add_rob_legend Logical indicating whether to add RoB legend. Default is FALSE.
#' @param exclude_high_rob Logical indicating whether to exclude high risk of bias studies 
#'   and refit the model. Default is FALSE.
#' @param font Character string specifying the font family to use throughout the plot.
#'   Default is NULL (uses system defaults). Common options include "Arial", "Times New Roman",
#'   "Helvetica", "Georgia", etc. The font must be available on your system.
#'
#' @return A patchwork object containing the complete forest plot with study information table,
#'   density plots, and effect size table.
#'
#' @details
#' The function performs several key operations:
#' \itemize{
#'   \item Validates input parameters and data structure
#'   \item Renames data columns based on the specified measure type
#'   \item Updates the brms model if necessary
#'   \item Handles Risk of Bias data if specified
#'   \item Creates either single or subgroup forest plots
#'   \item Combines tables and plots using patchwork
#' }
#'
#' For Risk of Bias functionality, the data must contain columns named "D1", "D2", "D3", and "Overall".
#' 
#' When \code{exclude_high_rob = TRUE}, studies with "High" overall risk of bias are removed
#' and the model is refitted on the filtered data.
#'
#' @examples
#' \dontrun{
#' # Basic usage with odds ratio
#' forest_plot <- bayes_forest(
#'   model = my_brms_model,
#'   data = my_data,
#'   measure = "OR",
#'   studyvar = author,
#'   year = year,
#'   c_n = control_n,
#'   i_n = intervention_n,
#'   c_event = control_events,
#'   i_event = intervention_events
#' )
#' 
#' # Subgroup analysis with mean difference
#' forest_plot_subgroup <- bayes_forest(
#'   model = my_brms_model,
#'   data = my_data,
#'   measure = "MD",
#'   studyvar = author,
#'   c_n = control_n,
#'   i_n = intervention_n,
#'   c_mean = control_mean,
#'   i_mean = intervention_mean,
#'   c_sd = control_sd,
#'   i_sd = intervention_sd,
#'   subgroup = TRUE,
#'   sort_subgroup_by = "effect"
#' )
#' }
#'
#' @seealso
#' \code{\link[brms]{brm}} for fitting Bayesian models,
#' \code{\link[metafor]{metafor}} for creating meta analysis models,
#' \code{\link[gt]{gt}} for making great tables,
#' \code{\link[patchwork]{patchwork}} for combining plots
#'
#' @author [Author Name]
#' @export
bayes_forest <- function(model,
                         data,
                         measure,
                         studyvar = NULL,
                         year = NULL,
                         c_n = NULL,
                         i_n = NULL,
                         c_event = NULL,
                         i_event = NULL,
                         c_mean = NULL,
                         i_mean = NULL,
                         c_sd = NULL,
                         i_sd = NULL,
                         c_time = NULL,
                         i_time = NULL,
                         sort_studies_by = "author",
                         subgroup = FALSE,
                         sort_subgroup_by = "alphabetical",
                         label_outcome = "Outcome",
                         label_control = "Control",
                         label_intervention = "Intervention",
                         title = NULL,
                         subtitle = NULL,
                         title_align = "left",
                         xlim = NULL,
                         x_breaks = NULL,
                         add_rope = FALSE,
                         rope_value = NULL,
                         rope_color = "grey50",
                         shrinkage_output = "density",
                         null_value = NULL,
                         color_palette = NULL,
                         color_study_posterior_null_left = "deepskyblue",
                         color_study_posterior_null_right = "violet",
                         color_study_posterior = "dodgerblue",
                         color_study_posterior_outline = "blue",
                         color_overall_posterior = "blue",
                         color_pointinterval = "purple",
                         color_shrinkage_outline = "purple",
                         color_shrinkage_fill = NULL,
                         plot_width = 4,
                         add_rob = FALSE,
                         rob_tool = c("rob2", "robins_i", "quadas2", "robins_e"),
                         add_rob_legend = FALSE,
                         exclude_high_rob = FALSE,
                         font = NULL) {
  
  # Input validation
  if (!inherits(model, "brmsfit")) {
    stop("Input must be a 'brmsfit' object")
  }
  
  if (is.null(data)) {
    stop("Data must be provided")
  }
  
  if (is.null(measure)) {
    stop("A measure must be entered: 'OR', 'HR', 'RR', 'IRR', 'MD', 'SMD'")
  }
  if (!measure %in% c("OR", "HR", "RR", "IRR", "MD", "SMD")) {
    stop("measure must be one of: 'OR', 'HR', 'RR', 'IRR', 'MD', 'SMD'")
  }
  
  if (!sort_studies_by %in% c("author", "year", "effect")) {
    stop("sort_studies_by must be one of 'author', 'year', or 'effect'")
  }
  
  # Get required variables for the measure
  if (measure %in% c("OR", "RR")) {
    data <- data |> 
      dplyr::rename(
        Author = {{studyvar}},
        Year = {{year}},
        N_Control = {{c_n}},
        N_Intervention = {{i_n}},
        Event_Control = {{c_event}},
        Event_Intervention = {{i_event}}
      )
    
  } else if (measure %in% c("MD", "SMD")) {
    data <- data |> 
      dplyr::rename(
        Author = {{studyvar}},
        Year = {{year}},
        N_Control = {{c_n}},
        N_Intervention = {{i_n}},
        Mean_Control = {{c_mean}},
        Mean_Intervention = {{i_mean}},
        SD_Control = {{c_sd}},
        SD_Intervention = {{i_sd}}
      )
    
  } else if (measure == "IRR") {
    data <- data |> 
      dplyr::rename(
        Author = {{studyvar}},
        Year = {{year}},
        Time_Control = {{c_time}},
        Time_Intervention = {{i_time}},
        Event_Control = {{c_event}},
        Event_Intervention = {{i_event}}
      )
  }
  
  # Update model
  model <- update_model(model, data, studyvar = Author)
  
  # Determine RoB variables in data
  rob_vars <- c("D1", "D2", "D3", "Overall")
  missing_vars <- setdiff(rob_vars, names(data))
  
  if ((length(missing_vars) > 0) && isTRUE(add_rob)) {
    stop("Risk of Bias columns must be provided for addition to the forest plot")
  }
  
  # Create columns to pass to internal functions so they won't throw an error
  if (length(missing_vars) > 0) {
    data <- data |> dplyr::mutate(
      D1 = NA_character_,
      Overall = NA_character_)
  }
  
  #if (!plot_output %in% c("density", "boxplot")) {
    #stop("shrinkage_output must be either 'density' or 'boxplot'")
  #}
  
  if (!shrinkage_output %in% c("density", "pointinterval")) {
    stop("shrinkage_output must be either 'density' or 'pointinterval'")
  }
  
  # Handle color_shrinkage_fill default
  if (is.null(color_shrinkage_fill)) {
    color_shrinkage_fill <- NA
  }
  
  # Handle exclude_high_rob if specified
  if (isTRUE(exclude_high_rob)) {
    if ("Overall" %in% names(data)) {
      data <- data |> dplyr::filter(Overall != "High" | is.na(Overall))
      
      # Refit the model on the filtered data
      message("Re-fitting model after excluding high risk of bias studies...")
      
      # Refit using the original formula and priors from the existing model
      model <- update(model, newdata = data, recompile = FALSE, refresh = 0)
    }
  }
  
  # Count distinct studies after RoB filtering
  num_studies <- data |> dplyr::distinct(Author) |> nrow()
  
  if (!subgroup) {
    if (num_studies < 2) {
      stop("Cannot run meta-analysis with fewer than 2 studies.")
    }
  } else {
    # Identify how many studies per subgroup
    subgroup_counts <- data |>
      dplyr::filter(!is.na(Subgroup)) |>
      dplyr::group_by(Subgroup) |>
      dplyr::summarise(n = dplyr::n_distinct(Author), .groups = "drop")
    
    # If no valid subgroups remain, warn user
    if (length(subgroup_counts) == 0) {
      warning("Subgroups with <2 studies will only show individual results and not pooled results for that subgroup. The study will still contribute to the overall effect.")
    }
  }
  
  # Create custom group order for subgroups if needed
  if (subgroup && "Subgroup" %in% names(data)) {
    if (is.character(sort_subgroup_by) && length(sort_subgroup_by) == 1) {
      if (sort_subgroup_by == "alphabetical") {
        subgroup_order <- c(sort(unique(data$Subgroup[!is.na(data$Subgroup)])), "Overall")
      } else if (sort_subgroup_by == "effect") {
        subgroup_effects <- data |>
          dplyr::group_by(Subgroup) |>
          dplyr::summarise(mean_effect = mean(yi, na.rm = TRUE), .groups = "drop") |>
          dplyr::arrange(mean_effect)
        subgroup_order <- c(subgroup_effects$Subgroup, "Overall")
      }
    } else if (is.character(sort_subgroup_by) && length(sort_subgroup_by) > 1) {
      # User supplied a custom group order directly
      subgroup_order <- c(sort_subgroup_by, "Overall")
    } else {
      stop("Invalid input to sort_subgroup_by. Must be 'alphabetical', 'effect', or a character vector of subgroup names.")
    }
  }
  
  # Handle different workflows for single vs subgroup plots
  if (!subgroup) {
    # Single Forest Plot Workflow
    # Step 1: Extract and organize posterior draws
    study.effect.draws <- forest.data_fn(
      data = data,
      model = model,
      subgroup = FALSE,
      sort_studies_by = sort_studies_by,
      subgroup_order = NULL
    )
    
    # Create the density plot
    study.plot <- study.density.plot_fn(
      df = study.effect.draws,
      model = model,
      measure = measure,
      subgroup = FALSE,
      subgroup_order = NULL,
      color_palette = color_palette,
      color_study_posterior = color_study_posterior,
      color_study_posterior_outline = color_study_posterior_outline,
      color_overall_posterior = color_overall_posterior,
      color_shrinkage_outline = color_shrinkage_outline,
      color_pointinterval = color_shrinkage_outline,
      color_shrinkage_fill = color_shrinkage_fill,
      label_control = label_control,
      label_intervention = label_intervention,
      shrinkage_output = shrinkage_output,
      xlim = xlim,
      x_breaks = x_breaks,
      null_value = null_value,
      add_rope = add_rope,
      rope_value = rope_value,
      rope_color = rope_color,
      font = font
    )
    
    # Step 3: Create summary data for tables
    forest.data.summary <- forest.data.summary_fn(
      spread_df = study.effect.draws,
      data = data,
      measure = measure,
      sort_studies_by = "author",
      subgroup = FALSE)
    
  } else {
    # Subgroup Forest Plot Workflow
    # Extract and organize posterior draws
    subgroup.effect.draws <- forest.data_fn(
      data = data,
      model = model,
      subgroup = TRUE,
      sort_studies_by = sort_studies_by,
      subgroup_order = subgroup_order
    )
    
    # Create subgroup summary
    forest.data.summary <- subgroup.effect.draws |>
      dplyr::mutate(Subgroup = factor(Subgroup, levels = subgroup_order)) |>
      dplyr::group_by(Subgroup) |>
      tidyr::nest() |>
      dplyr::rename(spread_df = data) |>
      dplyr::mutate(
        # Inject Subgroup back into each spread_df before passing to summary function
        spread_df = purrr::map2(spread_df, Subgroup, ~ dplyr::mutate(.x, Subgroup = .y)),
        subgroup.forest.summary = purrr::map(spread_df, ~ forest.data.summary_fn(
          spread_df = .x,
          data = data,  # full original data
          measure = measure,
          sort_studies_by = sort_studies_by))) |>
      tidyr::unnest(subgroup.forest.summary) |>
      dplyr::ungroup() |> 
      dplyr::select(-spread_df)
    
    if (isTRUE(subgroup) && measure %in% c("MD", "SMD")){
      forest.data.summary <- forest.data.summary |>
        dplyr::mutate(
          N_int = dplyr::case_when(
            Author == "Overall Effect" ~ as.character(sum(N_Intervention, na.rm = TRUE)),
            TRUE ~ N_int),
          N_ctrl = dplyr::case_when(
            Author == "Overall Effect" ~ as.character(sum(N_Control, na.rm = TRUE)),
            TRUE ~ N_ctrl),
          int_mean_sd = dplyr::case_when(
            Author == "Overall Effect" ~ NA,
            TRUE ~ int_mean_sd),
          ctrl_mean_sd = dplyr::case_when(
            Author == "Overall Effect" ~ NA,
            TRUE ~ ctrl_mean_sd))
    }
    
    if (isTRUE(subgroup) && measure %in% c("OR", "HR", "RR", "IRR")) {
      forest.data.summary <- forest.data.summary |>
        dplyr::mutate(
          Subgroup = dplyr::case_when(
            Author == "Overall Effect" ~ "Overall",
            TRUE ~ Subgroup),
          control_outcome_frac = dplyr::case_when(
            Author == "Overall Effect" & control_outcome_frac == "NA/NA" ~
              paste0(sum(Event_Control, na.rm = TRUE), "/", sum(N_Control, na.rm = TRUE)),
            TRUE ~ control_outcome_frac),
          int_outcome_frac = dplyr::case_when(
            Author == "Overall Effect" & int_outcome_frac == "NA/NA" ~
              paste0(sum(Event_Intervention, na.rm = TRUE), "/", sum(N_Intervention, na.rm = TRUE)),
            TRUE ~ int_outcome_frac))
    }
    
    
    if (isTRUE(exclude_high_rob)) {
      # If only 1 study, replace "Pooled Effect" with "No Pooled Effect"
      forest.data.summary <- forest.data.summary |>
        dplyr::left_join(subgroup_counts, by = "Subgroup") |>
        dplyr::mutate(
          Author = dplyr::case_when(
            n == 1 & Author %in% c("Pooled Effect", "Overall Effect") ~ "No Pooled Effect",
            TRUE ~ Author))
    } else {
      forest.data.summary
    }
    
    # Create plot data with subgroup (inserting spacers)
    subgroup.plot.data <- purrr::map_dfr(unique(subgroup.effect.draws$Subgroup), ~ {
      subgroup_data <- subgroup.effect.draws |> dplyr::filter(Subgroup == .x)
      spacer_row <- subgroup_data[1, ]
      spacer_row[1, ] <- NA
      spacer_row$Subgroup[1] <- .x
      spacer_row$Author[1] <- paste0("--- ", .x, " ---")
      dplyr::bind_rows(spacer_row, subgroup_data)
    })
    
    # Create lookup for author ordering
    seen_combos <- character(0)
    counter <- 1
    author_lookup <- purrr::map_dfr(1:nrow(subgroup.plot.data), ~ {
      current_combo <- paste(subgroup.plot.data$Subgroup[.x], subgroup.plot.data$Author[.x], sep = "_")
      
      if(!current_combo %in% seen_combos) {
        seen_combos <<- c(seen_combos, current_combo)  # Update the seen combos
        result <- data.frame(
          Subgroup = subgroup.plot.data$Subgroup[.x],
          Author = subgroup.plot.data$Author[.x],
          Author_ordered = counter
        )
        counter <<- counter + 1
        return(result)
      } else {
        return(NULL)
      }
    }) |>
      dplyr::filter(!is.null(Author_ordered))
    
    # Join back and create plot
    subgroup.plot.data <- subgroup.plot.data |>
      dplyr::left_join(author_lookup, by = c("Subgroup", "Author"))
    
    study.plot <- study.density.plot_fn(
      df = subgroup.plot.data,
      model = model,
      measure = measure,
      subgroup = subgroup,
      subgroup_order = subgroup_order,
      color_palette = color_palette,
      color_study_posterior = color_study_posterior,
      color_study_posterior_outline = color_study_posterior_outline,
      color_overall_posterior = color_overall_posterior,
      color_shrinkage_outline = color_shrinkage_outline,
      color_pointinterval = color_pointinterval,
      color_shrinkage_fill = color_shrinkage_fill,
      label_control = label_control,
      label_intervention = label_intervention,
      shrinkage_output = shrinkage_output,
      xlim = xlim,
      x_breaks = x_breaks,
      null_value = null_value,
      add_rope = add_rope,
      rope_value = rope_value,
      rope_color = rope_color,
      font = font
    )
  }
  
  # Create left table (study information)
  forest.table.left <- forest.table.left_fn(
    forest.data.summary = forest.data.summary,
    subgroup = subgroup,
    label_control = label_control,
    label_intervention = label_intervention,
    measure = measure,
    font = font
  )
  
  # Create right table (effect sizes and optionally RoB)
  forest.table.right <- forest.table.right_fn(
    df = forest.data.summary,
    subgroup = subgroup,
    add_rob = add_rob,
    measure = measure,
    font = font
  )
  
  # Combine everything using patchwork
  # Validate RoB legend parameters
  if (add_rob_legend == TRUE && add_rob == FALSE) {
    warning("Risk of Bias legend cannot be added when add_rob = FALSE. ",
            "Set add_rob = TRUE to include the RoB table and legend, ",
            "or set add_rob_legend = FALSE to suppress this warning.",
            call. = FALSE)
    add_rob_legend <- FALSE
  }
  
  # Combine everything using patchwork
  forest_plot <- patchwork_fn(
    table.left = forest.table.left,
    study.density.plot = study.plot,
    table.right = forest.table.right,
    plot_width = plot_width,
    title = title,
    subtitle = subtitle,
    title_align = title_align,
    add_rob_legend = add_rob_legend,
    rob_tool = rob_tool,
    font = font
  )
  
  return(forest_plot)
}
