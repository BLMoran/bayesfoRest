bayes_forest <- function(model,
                         sort_studies_by = "author",
                         subgroups = FALSE,
                         sort_subgroups_by = "alphabetical",
                         custom_group_order = NULL,
                         effectsize = NULL,
                         label_outcome = "Outcome",
                         label_control = "Control",
                         label_intervention = "Intervention",
                         xlim = NULL,
                         plot_output = "density",
                         shrinkage_output = "density",
                         density_fill_null = FALSE,
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
                         rob_tool = "rob2",
                         rob_legend = FALSE,
                         exclude_high_rob = FALSE,
                         data = NULL) {

  # Input validation
  if (!inherits(model, "brmsfit")) {
    stop("Input must be a 'brmsfit' object")
  }
  if (is.null(effectsize)) {
    stop("An effect size must be entered: 'OR', 'HR', 'MD', or 'SMD'")
  }
  if (is.null(data)) {
    stop("Data must be provided")
  }
  if (!sort_studies_by %in% c("author", "year", "effect")) {
    stop("sort_studies_by must be one of 'author', 'year', or 'effect'")
  }
  if (!effectsize %in% c("OR", "HR", "MD", "SMD")) {
    stop("effectsize must be one of 'OR', 'HR', 'MD', or 'SMD'")
  }
  if (!shrinkage_output %in% c("density", "pointinterval")) {
    stop("shrinkage_output must be either 'density' or 'pointinterval'")
  }

  # Handle color_shrinkage_fill default
  if (is.null(color_shrinkage_fill)) {
    color_shrinkage_fill <- NA
  }

  # Handle exclude_high_rob if specified
  if (exclude_high_rob && add_rob) {
    # Filter out studies with high risk of bias in Overall column
    if ("Overall" %in% names(data)) {
      data <- data |> dplyr::filter(Overall != "High" | is.na(Overall))
    }
  }

  # Create custom group order for subgroups if needed
  #custom_group_order <- NULL
  if (subgroups && "Subgroup" %in% names(data)) {
    if (sort_subgroups_by == "alphabetical") {
      custom_group_order <- c(sort(unique(data$Subgroup[!is.na(data$Subgroup)])), "Overall")
    } else if (sort_subgroups_by == "effect") {
      # Sort subgroups by their pooled effect size
      subgroup_effects <- data |>
        dplyr::group_by(Subgroup) |>
        dplyr::summarise(mean_effect = mean(yi, na.rm = TRUE), .groups = "drop") |>
        dplyr::arrange(mean_effect)
      custom_group_order <- c(subgroup_effects$Subgroup, "Overall")
    } else if (sort_subgroups_by == "custom") {
      custom_group_order <- c(custom_group_order, "Overall")
    }
  }

  # Handle different workflows for single vs subgroup plots
  if (!subgroups) {
    # Single Forest Plot Workflow
    # Step 1: Extract and organize posterior draws
    study.effect.draws <- forest.data_fn(
      data = data,
      model = model,
      subgroup = FALSE,
      sort_studies_by = sort_studies_by,
      custom_group_order = custom_group_order
    )

    # Step 2: Create the density plot
    study.plot <- study.density.plot_fn(
      df = study.effect.draws,
      effectsize = effectsize,
      sort_studies_by = sort_studies_by,
      subgroups = FALSE,
      custom_group_order = custom_group_order,
      color_study_posterior = color_study_posterior,
      color_study_posterior_outline = color_study_posterior_outline,
      color_overall_posterior = color_overall_posterior,
      color_shrinkage_outline = color_shrinkage_outline,
      color_pointinterval = color_shrinkage_outline,
      color_shrinkage_fill = color_shrinkage_fill,
      label_control = label_control,
      label_intervention = label_intervention,
      shrinkage_output = shrinkage_output
    )

    # Step 3: Create summary data for tables
    forest.data.summary <- forest.data.summary_fn(
      spread_df = study.effect.draws,
      data = data,
      sort_studies_by = sort_studies_by,
      #add_rob = add_rob,
      subgroups = FALSE
    )

  } else {
    # Subgroup Forest Plot Workflow
    # Step 1: Extract and organize posterior draws
    subgroup.effect.draws <- forest.data_fn(
      data = data,
      model = model,
      subgroup = TRUE,
      sort_studies_by = sort_studies_by,
      custom_group_order = custom_group_order
    )

    # Step 2: Create complex subgroup summary
    forest.data.summary <- subgroup.effect.draws |>
      dplyr::mutate(Subgroup = factor(Subgroup, levels = custom_group_order)) |>
      dplyr::group_by(Subgroup) |>
      tidyr::nest() |>
      dplyr::rename(spread_df = data) |>
      dplyr::mutate(
        # Inject Subgroup back into each spread_df before passing to summary function
        spread_df = purrr::map2(spread_df, Subgroup, ~ dplyr::mutate(.x, Subgroup = .y)),
        subgroup.forest.summary = purrr::map(spread_df, ~ forest.data.summary_fn(
          spread_df = .x,
          data = data,  # full original data
          #add_rob = add_rob,
          sort_studies_by = sort_studies_by
        ))
      ) |>
      tidyr::unnest(subgroup.forest.summary) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        Subgroup = dplyr::case_when(
          Author == "Overall Effect" ~ "Overall",
          TRUE ~ Subgroup),
        control_outcome_frac = dplyr::case_when(
          Author == "Overall Effect" & control_outcome_frac == "NA/NA" ~
            paste0(sum(Outcome_Control_Yes, na.rm = TRUE), "/", sum(N_Control, na.rm = TRUE)),
          TRUE ~ control_outcome_frac),
        int_outcome_frac = dplyr::case_when(
          Author == "Overall Effect" & int_outcome_frac == "NA/NA" ~
            paste0(sum(Outcome_Intervention_Yes, na.rm = TRUE), "/", sum(N_Intervention, na.rm = TRUE)),
          TRUE ~ int_outcome_frac)
      ) |>
      dplyr::select(-spread_df)

    # Step 3: Create plot data with subgroups (inserting spacers)
    subgroup.plot.data <- purrr::map_dfr(unique(subgroup.effect.draws$Subgroup), ~ {
      subgroup_data <- subgroup.effect.draws |> dplyr::filter(Subgroup == .x)
      spacer_row <- subgroup_data[1, ]
      spacer_row[1, ] <- NA
      spacer_row$Subgroup[1] <- .x
      spacer_row$Author[1] <- paste0("--- ", .x, " ---")
      dplyr::bind_rows(spacer_row, subgroup_data)
    })

    # Step 4: Create lookup for author ordering
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

    # Step 5: Join back and create plot
    subgroup.plot.data <- subgroup.plot.data |>
      dplyr::left_join(author_lookup, by = c("Subgroup", "Author"))

    study.plot <- study.density.plot_fn(
      df = subgroup.plot.data,
      effectsize = effectsize,
      sort_studies_by = sort_studies_by,
      subgroups = TRUE,
      custom_group_order = custom_group_order,
      color_study_posterior = color_study_posterior,
      color_study_posterior_outline = color_study_posterior_outline,
      color_overall_posterior = color_overall_posterior,
      color_shrinkage_outline = color_shrinkage_outline,
      color_pointinterval = color_pointinterval,
      color_shrinkage_fill = color_shrinkage_fill,
      label_control = label_control,
      label_intervention = label_intervention,
      shrinkage_output = shrinkage_output
    )
  }

  # Step 4: Create left table (study information)
  forest.table.left <- forest.table.left_fn(
    forest.data.summary = forest.data.summary,
    subgroups = subgroups,
    label_control = label_control,
    label_intervention = label_intervention,
    effectsize = effectsize
  )

  # Step 5: Create right table (effect sizes and optionally RoB)
  forest.table.right <- forest.table.right_fn(
    df = forest.data.summary,
    subgroups = subgroups,
    add_rob = add_rob,
    effectsize = effectsize
  )

  # Step 6: Combine everything using patchwork
  forest_plot <- patchwork_fn(
    table.left = forest.table.left,
    study.density.plot = study.plot,
    table.right = forest.table.right,
    plot_width = plot_width
  )

  # Return the complete forest plot
  return(forest_plot)
}

