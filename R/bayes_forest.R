#' Create Bayesian Forest Plot for Meta-Analysis
#'
#' @description
#' Creates a Bayesian forest plot for meta-analysis using a fitted Bayesian
#' model. Supports multiple backends (brms, MetaStan, bayesmeta, rstanarm,
#' rstan, cmdstanr, bayesma, RoBMA, R2jags) and both common-effect and
#' random-effect models.
#'
#' @param model A fitted Bayesian meta-analysis model object. Supported classes:
#'   \code{brmsfit} (brms), \code{meta_stan} (MetaStan), \code{bayesmeta},
#'   \code{stanreg} (rstanarm), \code{stanfit} (rstan), \code{CmdStanMCMC} (cmdstanr),
#'   \code{bayesma} (bayesma), \code{brma} (RoBMA), \code{rjags} (R2jags/JAGS).
#' @param data A data frame containing the study data used for the meta-analysis.
#' @param measure Character string specifying the effect measure. Must be one of:
#'   "OR" (Odds Ratio), "HR" (Hazard Ratio), "RR" (Risk Ratio),
#'   "IRR" (Incidence Rate Ratio), "MD" (Mean Difference), or "SMD" (Standardized Mean Difference).
#' @param studyvar Column name containing study identifiers/authors.
#' @param yearvar Column name containing publication years. Default is NULL.
#' @param n_ctrl Column name containing control group sample sizes.
#' @param n_int Column name containing intervention group sample sizes.
#' @param event_ctrl Column name containing control group event counts.
#' @param event_int Column name containing intervention group event counts.
#' @param mean_ctrl Column name containing control group means.
#' @param mean_int Column name containing intervention group means.
#' @param sd_ctrl Column name containing control group standard deviations.
#' @param sd_int Column name containing intervention group standard deviations.
#' @param time_ctrl Column name containing control group time periods.
#' @param time_int Column name containing intervention group time periods.
#' @param sort_studies_by Character string specifying how to sort studies.
#'   Options: "author" (default), "year", or "effect".
#' @param model_type Pooling type used to draw the plot: "random_effect" or
#'   "common_effect". If not supplied, it is detected from `model`.
#' @param re_min_k Optional minimum number of studies for a random-effect
#'   model. Between-study heterogeneity (tau) is poorly identified when there
#'   are few studies. If `re_min_k` is set and there are fewer studies than
#'   this, a random-effect `model_type` is switched to "common_effect" with a
#'   warning. `bayesma` models are refit as common-effect (overall and per
#'   subgroup); for other backends only the plot changes, as the pooled estimate
#'   still comes from the fitted random-effect model, so refit the model as
#'   common-effect yourself. Default is NULL (no minimum). An alternative is to
#'   keep `re_min_k = NULL` and fit the model with an informative prior on tau.
#' @param subgroup Logical indicating whether to create subgroup analysis. Default is FALSE.
#' @param subgroupvar Column name in `data` to use for subgroup analysis.
#' @param sort_subgroup_by Character string or vector specifying subgroup ordering.
#' @param label_outcome Character string for outcome label. Default is "Outcome".
#' @param label_control Character string for control group label. Default is "Control".
#' @param label_intervention Character string for intervention group label. Default is "Intervention".
#' @param title Character string for the plot title. Default is NULL (no title).
#' @param subtitle Character string for the plot subtitle. Default is NULL (no subtitle).
#' @param title_align Character string specifying title alignment.
#' @param xlim Numeric vector of length 2 specifying x-axis limits. Default is NULL.
#' @param x_breaks Numeric vector specifying custom x-axis break points. Default is NULL.
#' @param x_label Character string for the x-axis label. Default is NULL
#'   (measure-specific label).
#' @param effect_label Character string used for the effect columns in the
#'   right-hand table. Default is NULL (uses `measure`).
#' @param plot_output Character string specifying how study and pooled effects
#'   are drawn: "density" (default) or "boxplot". Boxplots are horizontal: the
#'   box spans the interquartile range, the line marks the median, and the
#'   whiskers span the 95% interval. With "boxplot", `shrinkage_output` must be
#'   "pointinterval" (it is set automatically if not supplied).
#' @param shrinkage_output Character string specifying shrinkage visualization.
#'   Options: "density" (default) or "pointinterval". Ignored for common-effect models.
#' @param null_value Numeric value specifying x-axis value of the null line. Default is NULL.
#' @param add_rope Logical indicating whether to add ROPE. Default is FALSE.
#' @param rope_value Numeric vector of length 2 specifying ROPE range.
#' @param rope_color Color for ROPE shading.
#' @param color_palette Character vector of colors for the plot. Default is NULL.
#' @param color_study_posterior_null_left Color for left side of null study posteriors.
#' @param color_study_posterior_null_right Color for right side of null study posteriors.
#' @param color_study_posterior Color for study posterior densities.
#' @param color_study_posterior_outline Color for study posterior outlines.
#' @param color_overall_posterior Color for overall posterior.
#' @param color_pointinterval Color for point intervals.
#' @param color_shrinkage_outline Color for shrinkage plot outlines.
#' @param color_shrinkage_fill Color for shrinkage plot fill.
#' @param split_color_by_null Logical. Split posterior densities by null value.
#'   For boxplots, each box is coloured by which side of the null its median falls.
#' @param color_favours_control Colour for density regions favouring the control group.
#' @param color_favours_intervention Colour for density regions favouring the intervention group.
#' @param add_arm_labels Logical. Whether to show "Favours Control/Intervention" labels. Default is TRUE.
#' @param reverse_arm_labels Logical. Whether to swap the arm label positions. Default is FALSE.
#' @param add_pred Logical. Whether to add a predictive distribution for a new study. Default is FALSE.
#'   Only applicable to random-effect models.
#' @param add_pred_subgroup Logical. Whether to add prediction to each subgroup. Default is FALSE.
#' @param pred_output Character string specifying prediction visualization.
#'   Options: "density" (default) or "pointinterval".
#' @param color_pred_posterior Color for prediction density fill. Default is "forestgreen".
#' @param color_pred_outline Color for prediction density outline. Default is "darkgreen".
#' @param color_pred_pointinterval Color for prediction point interval. Default is "forestgreen".
#' @param plot_width Numeric value specifying the relative width of the plot component. Default is 4.
#' @param add_rob Logical indicating whether to add Risk of Bias assessment. Default is FALSE.
#' @param rob_tool Risk of bias tool: "rob2" (default), "robins_i", "robins_e",
#'   "quadas2", "clarity", "rob_nma", "probast", or "probast_ai". See
#'   [rob_plot()] for the domain columns each tool expects.
#' @param add_rob_legend Logical indicating whether to add RoB legend. Default is FALSE.
#' @param exclude_high_rob Logical indicating whether to exclude studies with a
#'   high overall risk of bias ("High", "Serious", "Critical", "Very high" or
#'   "Definitely high") and refit the model. Requires an overall RoB column.
#' @param font Character string specifying the font family. Default is NULL.
#'
#' @return A patchwork object containing the complete forest plot.
#'
#' @export
bayes_forest <- function(model,
                         data,
                         measure,
                         studyvar = NULL,
                         yearvar = NULL,
                         n_ctrl = NULL,
                         n_int = NULL,
                         event_ctrl = NULL,
                         event_int = NULL,
                         mean_ctrl = NULL,
                         mean_int = NULL,
                         sd_ctrl = NULL,
                         sd_int = NULL,
                         time_ctrl = NULL,
                         time_int = NULL,
                         sort_studies_by = "author",
                         model_type = c("random_effect", "common_effect"),
                         re_min_k = NULL,
                         subgroup = FALSE,
                         subgroupvar = NULL,
                         sort_subgroup_by = "alphabetical",
                         label_outcome = "Outcome",
                         label_control = "Control",
                         label_intervention = "Intervention",
                         title = NULL,
                         subtitle = NULL,
                         title_align = "left",
                         xlim = NULL,
                         x_breaks = NULL,
                         x_label = NULL,
                         effect_label = NULL,
                         add_rope = FALSE,
                         rope_value = NULL,
                         rope_color = "grey50",
                         plot_output = c("density", "boxplot"),
                         shrinkage_output = c("density", "pointinterval"),
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
                         split_color_by_null = FALSE,
                         color_favours_control = "firebrick",
                         color_favours_intervention = "dodgerblue",
                         add_arm_labels = TRUE,
                         reverse_arm_labels = FALSE,
                         add_pred = FALSE,
                         add_pred_subgroup = FALSE,
                         pred_output = c("density", "pointinterval"),
                         color_pred_posterior = "forestgreen",
                         color_pred_outline = "darkgreen",
                         color_pred_pointinterval = "forestgreen",
                         plot_width = 4,
                         add_rob = FALSE,
                         rob_tool = "rob2",
                         add_rob_legend = FALSE,
                         exclude_high_rob = FALSE,
                         font = NULL) {

  # Input validation
  validate_inputs(
    model = model,
    data = data,
    measure = measure,
    studyvar = {{studyvar}},
    yearvar = {{yearvar}},
    subgroup = subgroup,
    subgroupvar = {{subgroupvar}},
    sort_studies_by = sort_studies_by
  )

  # Boxplots only pair with point-interval shrinkage estimates
  plot_output <- rlang::arg_match(plot_output)
  if (plot_output == "boxplot") {
    if (!missing(shrinkage_output) && identical(shrinkage_output, "density")) {
      cli::cli_abort(c(
        "{.arg shrinkage_output} must be {.val pointinterval} when {.arg plot_output} = {.val boxplot}.",
        "i" = "Drop {.arg shrinkage_output} or set it to {.val pointinterval}."
      ))
    }
    shrinkage_output <- "pointinterval"
  } else {
    shrinkage_output <- rlang::arg_match(shrinkage_output)
  }

  # Auto-detect model type (or use user override)
  model_info <- detect_model_info(model)
  model_type <- if (missing(model_type)) model_info$pooling else rlang::arg_match(model_type)

  # Rename data columns based on measure type
  if (measure %in% c("OR", "RR")) {
    data <- data |>
      dplyr::rename(
        Author = {{studyvar}},
        Year = {{yearvar}},
        N_Control = {{n_ctrl}},
        N_Intervention = {{n_int}},
        Event_Control = {{event_ctrl}},
        Event_Intervention = {{event_int}}
      )

  } else if (measure %in% c("MD", "SMD")) {
    have_effect <- all(c("yi", "vi") %in% names(data))
    data <- data |>
      dplyr::rename(
        Author = {{studyvar}},
        Year = {{yearvar}},
        N_Control = {{n_ctrl}},
        N_Intervention = {{n_int}})
    if (!have_effect) {
      # only need arm mean/SD when yi/vi must be computed
      data <- data |>
        dplyr::rename(
          Mean_Control = {{mean_ctrl}},
          Mean_Intervention = {{mean_int}},
          SD_Control = {{sd_ctrl}},
          SD_Intervention = {{sd_int}})
    }
  } else if (measure == "IRR") {
    data <- data |>
      dplyr::rename(
        Author = {{studyvar}},
        Year = {{yearvar}},
        Time_Control = {{time_ctrl}},
        Time_Intervention = {{time_int}},
        Event_Control = {{event_ctrl}},
        Event_Intervention = {{event_int}}
      )
  }

  # Handle subgroup variable renaming
  if (isTRUE(subgroup)) {
    if (!rlang::quo_is_null(rlang::enquo(subgroupvar))) {
      subgroup_var_name <- rlang::as_name(rlang::ensym(subgroupvar))
      if (subgroup_var_name != "Subgroup") {
        data <- data |>
          dplyr::rename(Subgroup = {{subgroupvar}})
      }
    }
  }

  # Disambiguate duplicate Author names
  data <- make_authors_unique(data)

  # Compute yi/vi if not present (needed for one-stage models)
  data <- compute_yi_vi(data, measure)

  # Risk of bias: validate and standardise judgements for the chosen tool
  rob_spec_obj <- NULL
  rob_cols <- character()
  if (isTRUE(add_rob) || isTRUE(exclude_high_rob)) {
    rob_spec_obj <- as_rob_spec(rob_tool)
    data <- prepare_rob_data(data, rob_spec_obj)
    rob_cols <- rob_columns(data, rob_spec_obj)
  }

  # Handle color_shrinkage_fill default
  if (is.null(color_shrinkage_fill)) {
    color_shrinkage_fill <- NA
  }

  # Handle exclude_high_rob if specified
  if (isTRUE(exclude_high_rob)) {
    overall_col <- dplyr::intersect(names(rob_spec_obj@overall), names(data))[1]
    if (is.na(overall_col)) {
      cli::cli_abort("{.arg exclude_high_rob} requires an overall risk of bias column; {rob_spec_obj@label} has none.")
    }
    data <- data |> dplyr::filter(!.data[[overall_col]] %in% rob_high_levels())

    message("Re-fitting model after excluding high risk of bias studies...")
    refit <- refit_for_subgroup(model, data)
    if (!is.null(refit)) {
      model <- refit
    } else {
      rlang::warn("Could not refit model after RoB exclusion. Using original model.")
    }
  }

  # Match pred_output argument
  pred_output <- match.arg(pred_output)

  # Count distinct studies after RoB filtering
  num_studies <- data |> dplyr::distinct(Author) |> nrow()

  if (!subgroup) {
    if (num_studies < 2) {
      stop("Cannot run meta-analysis with fewer than 2 studies.")
    }
  } else {
    subgroup_counts <- data |>
      dplyr::filter(!is.na(Subgroup)) |>
      dplyr::group_by(Subgroup) |>
      dplyr::summarise(n = dplyr::n_distinct(Author), .groups = "drop")

    if (length(subgroup_counts) == 0) {
      warning("Subgroups with <2 studies will only show individual results and not pooled results for that subgroup. The study will still contribute to the overall effect.")
    }
  }

  # Enforce a minimum number of studies for random-effect models
  requested_type <- model_type
  model_type <- enforce_re_min_k(model_type, re_min_k, data, "Author")
  if (requested_type == "random_effect" && model_type == "common_effect" &&
      model_info$pooling == "random_effect") {
    if (inherits(model, "bayesma")) {
      refit <- refit_for_subgroup(model, data, re_min_k = re_min_k)
      if (!is.null(refit)) {
        model <- refit
      } else {
        rlang::warn("Could not refit bayesma model as common-effect. Using original model.")
      }
    } else {
      cli::cli_warn(c(
        "!" = "The fitted {.cls {class(model)[1]}} model is random-effect; only the plot is drawn as common-effect.",
        "i" = "Refit the model as common-effect for a common-effect pooled estimate,
               or keep {.arg re_min_k} = NULL and use an informative prior on tau."
      ))
    }
  }

  if (!is.null(re_min_k) && isTRUE(subgroup) && model_type == "random_effect" &&
      !inherits(model, "bayesma")) {
    small_subgroups <- subgroup_counts |>
      dplyr::filter(n < re_min_k) |>
      dplyr::pull(Subgroup)
    if (length(small_subgroups) > 0) {
      cli::cli_warn(c(
        "!" = paste0(
          "Subgroup{?s} {.val {small_subgroups}} ha{?s/ve} fewer than ",
          re_min_k, " studies; tau is poorly identified."),
        "i" = "Subgroup pooled effects are still random-effect. Consider an informative prior on tau."
      ))
    }
  }

  # Extract pooled effect summary once (used by density plot for reference lines)
  pooled_effect <- extract_pooled_summary(model)

  # Create subgroup order if needed
  subgroup_order <- if (isTRUE(subgroup)) {
    get_subgroup_order(data, sort_subgroup_by)
  } else {
    NULL
  }

  # Handle different workflows for single vs subgroup plots
  if (!subgroup) {
    # Single Forest Plot Workflow
    study.effect.draws <- forest.data_fn(
      data = data,
      model = model,
      subgroup = FALSE,
      sort_studies_by = sort_studies_by,
      subgroup_order = NULL,
      effect_type = model_type,
      add_pred = add_pred
    )

    study.plot <- study.density.plot_fn(
      df = study.effect.draws,
      pooled_effect = pooled_effect,
      measure = measure,
      effect_type = model_type,
      subgroup = FALSE,
      subgroup_order = NULL,
      color_palette = color_palette,
      color_study_posterior = color_study_posterior,
      color_study_posterior_outline = color_study_posterior_outline,
      color_overall_posterior = color_overall_posterior,
      color_shrinkage_outline = color_shrinkage_outline,
      color_pointinterval = color_shrinkage_outline,
      color_shrinkage_fill = color_shrinkage_fill,
      split_color_by_null = split_color_by_null,
      color_favours_control = color_favours_control,
      color_favours_intervention = color_favours_intervention,
      label_control = label_control,
      label_intervention = label_intervention,
      add_arm_labels = add_arm_labels,
      reverse_arm_labels = reverse_arm_labels,
      plot_output = plot_output,
      shrinkage_output = shrinkage_output,
      add_pred = add_pred,
      pred_output = pred_output,
      color_pred_posterior = color_pred_posterior,
      color_pred_outline = color_pred_outline,
      color_pred_pointinterval = color_pred_pointinterval,
      xlim = xlim,
      x_breaks = x_breaks,
      x_label = x_label,
      null_value = null_value,
      add_rope = add_rope,
      rope_value = rope_value,
      rope_color = rope_color,
      font = font
    )

    forest.data.summary <- forest.data.summary_fn(
      spread_df = study.effect.draws,
      data = data,
      measure = measure,
      sort_studies_by = sort_studies_by,
      subgroup = FALSE,
      effect_type = model_type,
      rob_cols = rob_cols)

  } else {
    # Subgroup Forest Plot Workflow
    subgroup.effect.draws <- forest.data_fn(
      data = data,
      model = model,
      subgroup = TRUE,
      sort_studies_by = sort_studies_by,
      subgroup_order = subgroup_order,
      effect_type = model_type,
      add_pred = add_pred,
      add_pred_subgroup = add_pred_subgroup,
      re_min_k = re_min_k
    )

    # Create subgroup summary
    forest.data.summary <- subgroup.effect.draws |>
      dplyr::mutate(Subgroup = factor(Subgroup, levels = subgroup_order)) |>
      dplyr::group_by(Subgroup) |>
      tidyr::nest() |>
      dplyr::rename(spread_df = data) |>
      dplyr::mutate(
        spread_df = purrr::map2(spread_df, Subgroup, ~ dplyr::mutate(.x, Subgroup = .y)),
        subgroup.forest.summary = purrr::map(spread_df, ~ forest.data.summary_fn(
          spread_df = .x,
          data = data,
          measure = measure,
          sort_studies_by = sort_studies_by,
          effect_type = model_type,
          rob_cols = rob_cols))) |>
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
        seen_combos <<- c(seen_combos, current_combo)
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
      pooled_effect = pooled_effect,
      measure = measure,
      effect_type = model_type,
      subgroup = subgroup,
      subgroup_order = subgroup_order,
      color_palette = color_palette,
      color_study_posterior = color_study_posterior,
      color_study_posterior_outline = color_study_posterior_outline,
      color_overall_posterior = color_overall_posterior,
      color_shrinkage_outline = color_shrinkage_outline,
      color_pointinterval = color_pointinterval,
      color_shrinkage_fill = color_shrinkage_fill,
      split_color_by_null = split_color_by_null,
      color_favours_control = color_favours_control,
      color_favours_intervention = color_favours_intervention,
      label_control = label_control,
      label_intervention = label_intervention,
      add_arm_labels = add_arm_labels,
      reverse_arm_labels = reverse_arm_labels,
      plot_output = plot_output,
      shrinkage_output = shrinkage_output,
      add_pred = add_pred,
      #add_pred_subgroup = add_pred_subgroup,
      pred_output = pred_output,
      color_pred_posterior = color_pred_posterior,
      color_pred_outline = color_pred_outline,
      color_pred_pointinterval = color_pred_pointinterval,
      xlim = xlim,
      x_breaks = x_breaks,
      x_label = x_label,
      null_value = null_value,
      add_rope = add_rope,
      rope_value = rope_value,
      rope_color = rope_color,
      font = font
    )
  }

  # Filter out Prediction rows from summary if they shouldn't be there
  # (density plot handles its own filtering via add_pred flag)

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
    rob_cols = rob_cols,
    spec = rob_spec_obj,
    measure = measure,
    font = font,
    effect_label = effect_label
  )

  # Validate RoB legend parameters
  if (add_rob_legend == TRUE && add_rob == FALSE) {
    warning("Risk of Bias legend cannot be added when add_rob = FALSE. ",
            "Set add_rob = TRUE to include the RoB table and legend, ",
            "or set add_rob_legend = FALSE to suppress this warning.",
            call. = FALSE)
    add_rob_legend <- FALSE
  }

  rob_legend <- if (isTRUE(add_rob_legend)) {
    levels_present <- data |>
      dplyr::select(dplyr::all_of(rob_cols)) |>
      purrr::map(as.character) |>
      purrr::list_c() |>
      unique()
    rob_legend_fn(
      rob_spec_obj,
      domains_present = dplyr::intersect(rob_cols, names(rob_spec_obj@domains)),
      levels_present = levels_present,
      font = font)
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
    rob_legend = rob_legend,
    font = font
  )

  return(forest_plot)
}
