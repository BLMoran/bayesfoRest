#' Internal function to sort studies
#'
#' @noRd
sort.studies.fn <- function(.data, sort_studies_by = NULL) {
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


#' Internal function to get certain properties dependent on measure
#'
#' @noRd
get_measure_properties <- function(measure) {
  switch(measure,
         "OR" = list(
           null_value = 1,
           log_scale = TRUE,
           #x_breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 7),
           x_label = "Odds Ratio (log scale)",
           data_cols = c("Event_Control", "N_Control", "Event_Intervention", "N_Intervention")
         ),
         "HR" = list(
           null_value = 1,
           log_scale = TRUE,
           #x_breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 7),
           x_label = "Hazard Ratio (log scale)",
           data_cols = c("Event_Control", "N_Control", "Event_Intervention", "N_Intervention")
         ),
         "RR" = list(
           null_value = 1,
           log_scale = TRUE,
           #x_breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 7),
           x_label = "Risk Ratio (log scale)",
           data_cols = c("Event_Control", "N_Control", "Event_Intervention", "N_Intervention")
         ),
         "IRR" = list(
           null_value = 1,
           log_scale = TRUE,
           #x_breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 7),
           x_label = "Incident Rate Ratio (log scale)",
           data_cols = c("Event_Control", "N_Control", "Time_Control", "Event_Intervention", "N_Intervention", "Time_Intervention")
         ),
         "MD" = list(
           null_value = 0,
           log_scale = FALSE,
           #x_breaks = NULL,
           x_label = "Mean Difference",
           data_cols = c("Mean_Control", "SD_Control", "N_Control", "Mean_Intervention", "SD_Intervention", "N_Intervention")
         ),
         "SMD" = list(
           null_value = 0,
           log_scale = FALSE,
           #x_breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2),
           x_label = "Standardised Mean Difference",
           data_cols = c("Mean_Control", "SD_Control", "N_Control", "Mean_Intervention", "SD_Intervention", "N_Intervention")
         ),
         stop("Effect size must be one of: 'OR', 'HR', 'RR', 'IRR', 'MD', 'SMD'")
  )
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
  model <- update(
    object = model,
    formula = new_formula,
    newdata = data,
    recompile = FALSE
  )
  return(model)
}



#' Internal function to extract draws from the posterior
#'
#' @noRd
forest.data_fn <- function(data,
                           model,
                           subgroup = FALSE,
                           sort_studies_by = "author",
                           subgroup_order = NULL) {
  
  if (subgroup == FALSE && "Subgroup" %in% names(data)) {
    data <- data |> dplyr::select(-Subgroup)
  }
  if (subgroup == FALSE) {
    # No subgroup, evaluate as a whole
    study.draws <- tidybayes::spread_draws(model, r_Author[Author, ], b_Intercept) |>
      dplyr::mutate(b_Intercept = r_Author + b_Intercept)
    pooled.draws <- tidybayes::spread_draws(model, b_Intercept, sd_Author__Intercept) |>
      dplyr::mutate(Author = "Pooled Effect")
    effect.draws <- dplyr::bind_rows(study.draws, pooled.draws) |>
      dplyr::mutate(Author = stringr::str_replace_all(Author, "\\.", " ")) |> 
      dplyr::ungroup() |>
      dplyr::left_join(dplyr::select(data, Author, Year, yi, vi), by = "Author") |>
      sort.studies.fn(sort_studies_by)
  } else {
    # With subgroup
    subgroup_df <- data |>
      dplyr::group_by(Subgroup) |>
      tidyr::nest() |>
      dplyr::mutate(
        subgroup_model = purrr::map(data, ~ update(model, newdata = .x)),
        # Calculate study count for each subgroup
        study_count = purrr::map_int(data, nrow)
      )
    
    # Create overall effect draws with matching structure to subgroup effect.draws
    overall.effect.draws <- tidybayes::spread_draws(model, b_Intercept, sd_Author__Intercept) |>
      dplyr::mutate(
        Author = "Overall Effect",
        Subgroup = "Overall",
        r_Author = 0,           # No random effect for overall
        Year = NA_integer_,     # No specific year
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
          
          combined <- dplyr::bind_rows(study, pooled) |>
            dplyr::mutate(Author = stringr::str_replace_all(Author, "\\.", " ")) |> 
            dplyr::ungroup() |>
            dplyr::left_join(dplyr::select(..2, Author, Year, yi, vi), by = "Author") |>
            sort.studies.fn(sort_studies_by)
          return(combined)
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
        dplyr::mutate(Subgroup =  dplyr::case_when(
          is.na(Subgroup) & Author == "Overall Effect" ~ "Overall",
          .default = Subgroup
        ))
    }
  }
  return(effect.draws)
}

#' Internal function to generate plot section of forest plot
#'
#' @noRd
study.density.plot_fn  <- function(df,
                                   model,
                                   measure = measure,
                                   subgroup = FALSE,
                                   subgroup_order = NULL,
                                   color_palette = NULL,
                                   color_study_posterior = "dodgerblue",
                                   color_study_posterior_outline = "blue",
                                   color_pooled_posterior = "blue",
                                   color_overall_posterior = "darkblue",
                                   color_shrinkage_outline = "purple",
                                   color_pointinterval = "purple",
                                   color_shrinkage_fill = NA,
                                   label_control = "Control",
                                   label_intervention = "Intervention",
                                   shrinkage_output = "density",
                                   xlim = NULL,
                                   x_breaks = NULL,
                                   null_value = NULL,
                                   add_rope = FALSE,
                                   rope_value = NULL,
                                   rope_color = "grey50",
                                   font = NULL){
  
  # Filter out "No Pooled Effect" rows at the beginning
  df <- df |> dplyr::filter(Author != "No Pooled Effect")
  
  # Apply color palette if provided
  if (!is.null(color_palette)) {
    # Extract colors from the specified palette
    # Assuming we need at least 6 colors for the various elements
    palette_colors <- paletteer::paletteer_d(color_palette, n = 6)
    
    # Override default colors with palette colors
    color_study_posterior <- palette_colors[2]
    color_study_posterior_outline <- palette_colors[1]
    color_pooled_posterior <- palette_colors[3]
    color_overall_posterior <- palette_colors[5]
    color_shrinkage_outline <- palette_colors[6]
    color_pointinterval <- palette_colors[6]
    
    # If you want to use a fill color from the palette instead of NA
    if (!is.na(color_shrinkage_fill)) {
      color_shrinkage_fill <- scales::alpha(palette_colors[6], 0.3)
    }
  }
  
  # Get effect size properties
  props <- get_measure_properties(measure) 
  
  # Determine x_breaks to use: custom if provided, otherwise use measure defaults
  breaks <- if (!is.null(x_breaks)) x_breaks else waiver()
  
  # Determine null line to use: custom if provided, otherwise use measure default
  null_value <- if (!is.null(null_value)) null_value else props$null_value
  
  # Process ROPE values if provided
  if (isTRUE(add_rope)) {
    if (is.null(rope_value)) {
      # Use Kruschke's recommended defaults if no rope_value provided
      if (measure %in% c("OR", "HR", "RR", "IRR")) {
        # For ratio measures: ROPE around 1.0 with ±10% range (0.9 to 1.1)
        rope_range <- c(0.9, 1.1)
      } else if (measure == "SMD") {
        # For standardized mean difference: Cohen's small effect size ±0.1
        rope_range <- c(-0.1, 0.1)
      } else if (measure == "MD") {
        # For mean difference: no universal default, require user specification
        stop("For mean differences, rope_value must be specified as it depends on the measurement scale")
      }
    } else {
      if (length(rope_value) == 1) {
        # Single value: create symmetric range around null line
        rope_range <- c(null_line_to_use - rope_value, null_line_to_use + rope_value)
      } else if (length(rope_value) == 2) {
        # Two values: use as explicit range
        rope_range <- sort(rope_value)  # Ensure correct order
      } else {
        stop("rope_value must be either a single value or a vector of length 2")
      }
    }
  }
  
  # Optimise data structure: separate study-level data from posterior draws
  if (isTRUE(subgroup)){
    study.effects <- df |> dplyr::distinct(Author_ordered, yi, vi) |>
      dplyr::mutate(
        Author_ordered = factor(Author_ordered, levels = as.character(1:max(Author_ordered, na.rm = TRUE))),
        Author = forcats::fct_rev(Author_ordered))
    
    posterior.draws <- df |>
      dplyr::mutate(
        Author_pooled = Author,
        Author_ordered = factor(Author_ordered, levels = as.character(1:max(Author_ordered, na.rm = TRUE))),
        Author = forcats::fct_rev(Author_ordered))
    
  } else if (isFALSE(subgroup)){
    study.effects <- df |> dplyr::distinct(Author, yi, vi)
    posterior.draws <- df
  }
  
  if (isTRUE(props$log_scale)){
    x.min <- 0.1
    x.max <- ceiling(max(exp(df$yi + 1.96 * sqrt(df$vi)), na.rm = TRUE))
  } else if (isFALSE(props$log_scale)){
    x.min <- floor(min(df$yi - 1.96 * sqrt(df$vi), na.rm = TRUE))
    x.max <- ceiling(max(df$yi + 1.96 * sqrt(df$vi), na.rm = TRUE))
  }
  
  # Set xlim to given value or calculated range
  if (!is.null(xlim)){
    calc_xlim <- xlim
  } else {
    calc_xlim <- c(x.min, x.max)
  }
  
  if (isTRUE(props$log_scale)){
    posterior.draws <- posterior.draws |> 
      dplyr::mutate(
        x_studies = exp(b_Intercept))
    
    study.effects <- study.effects |> 
      dplyr::mutate(
        xdist = exp(distributional::dist_normal(mean = yi, sd = sqrt(vi))))
  } else {
    posterior.draws <- posterior.draws |> 
      dplyr::mutate(
        x_studies = b_Intercept)
    
    study.effects <- study.effects |> 
      dplyr::mutate(
        xdist = distributional::dist_normal(mean = yi, sd = sqrt(vi)))
  }
  
  study.density.plot <-
    ggplot2::ggplot(ggplot2::aes(y = Author), data = posterior.draws) +
    # Add ROPE shading first (so it appears behind other elements)
    {if (isTRUE(add_rope)) {
      ggplot2::annotate("rect", 
                        xmin = rope_range[1], xmax = rope_range[2], 
                        ymin = -Inf, ymax = Inf, 
                        fill = scales::alpha(rope_color, 0.3), 
                        color = NA)
    }} +
    ggdist::stat_slab(ggplot2::aes(xdist = xdist),
                      slab_linewidth = 0.5, alpha = 0.7, limits = calc_xlim, height = 0.9,
                      data = study.effects, colour = color_study_posterior_outline, fill = color_study_posterior) +
    ggdist::stat_slab(ggplot2::aes(x = x_studies, y = Author),
                      data = posterior.draws |> dplyr::filter(if (isTRUE(subgroup)) {Author_pooled == "Pooled Effect"}
                                                              else {Author == "Pooled Effect"}),
                      fill = color_pooled_posterior, height = 0.9, normalize = "panels") +
    # Add overall effect slab when subgroup = TRUE
    {if (isTRUE(subgroup)) {
      ggdist::stat_slab(ggplot2::aes(x = x_studies, y = Author),
                        data = posterior.draws |> dplyr::filter(Author_pooled == "Overall Effect"),
                        fill = color_overall_posterior, height = 0.9, normalize = "panels")
    }} +
    ggplot2::guides(alpha = "none", fill = "none")+
    # Add null value
    ggplot2::geom_vline(xintercept = null_value, color = "black", linewidth = 1) +
    ggplot2::coord_cartesian(xlim= calc_xlim, clip = "off") +
    ggplot2::theme_light() +
    ggplot2::theme(
      axis.ticks.y = ggplot2::element_blank(), 
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(vjust=-0.5, family = font), 
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0,0,0,0), 
      panel.border = ggplot2::element_blank(),
      axis.line.x.top = ggplot2::element_line(color = "grey60", linewidth = 0.75),
      axis.line.x.bottom = ggplot2::element_line(color = "black", linewidth = 0.75),
      axis.text.x.top = ggplot2::element_blank(), 
      axis.ticks.x.top = ggplot2::element_blank(),
      axis.text.x.bottom = ggplot2::element_text(colour = "black", family = font)) +
    ggplot2::guides(x.sec = "axis", y.sec = "axis") +
    ggplot2::annotation_custom(grid::textGrob(
      label = paste(" Favours\n", label_intervention),
      x = grid::unit(-0.01, "npc"), y = grid::unit(1.02, "npc") , just = c("left", "bottom"),
      gp = grid::gpar(col = "grey30", fontsize = 10, fontfamily = font)),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
    ggplot2::annotation_custom(grid::textGrob(
      label = paste(" Favours\n", label_control),
      x = grid::unit(1, "npc"), y = grid::unit(1.02, "npc") , just = c("right", "bottom"),
      gp = grid::gpar(col = "grey30", fontsize = 10, fontfamily = font)),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
    ggplot2::labs(x= props$x_label) +
    ggplot2::ylab(NULL) +
    ggplot2::geom_hline(yintercept = 1, color = "black", linewidth = 0.75)
  
  if (isTRUE(props$log_scale)){
    study.density.plot <- study.density.plot +
      ggplot2::scale_x_log10(breaks = breaks, expand = c(0, 0), limits = calc_xlim)+
      ggplot2::geom_vline(xintercept = exp(brms::fixef(model)[1, 1]), color = "grey60", linewidth = 1) +
      ggplot2::geom_vline(xintercept = exp(brms::fixef(model)[1, 3:4]), color = "grey60", linetype = 2)
  } else {
    study.density.plot <- study.density.plot +
      ggplot2::scale_x_continuous(breaks = breaks, expand = c(0, 0), limits = calc_xlim) +
      ggplot2::geom_vline(xintercept = brms::fixef(model)[1, 1], color = "grey60", linewidth = 1) +
      ggplot2::geom_vline(xintercept = brms::fixef(model)[1, 3:4], color = "grey60", linetype = 2)
  }
  
  if(isFALSE(subgroup)){
    study.density.plot <- study.density.plot + ggplot2::scale_y_discrete(expand = c(0, 0), limits = rev)
  } else if (isTRUE(subgroup)){
    study.density.plot <- study.density.plot + ggplot2::scale_y_discrete(expand = c(0, 0))
  }
  
  if(is.null(shrinkage_output))"density" else shrinkage_output
  
  if (shrinkage_output == "density") {
    study.density.plot <- study.density.plot +
      ggdist::stat_slab(
        ggplot2::aes(x = x_studies, y = Author),
        data = posterior.draws |> dplyr::filter(if (isTRUE(subgroup)) {Author_pooled != "Pooled Effect" & Author_pooled != "Overall Effect"}
                                                else {Author != "Pooled Effect"}),
        linewidth = 0.5,
        scale = 0.6,
        normalize = "panels",
        color = color_shrinkage_outline,
        fill = color_shrinkage_fill,
        limits = calc_xlim
      )
  } else if (shrinkage_output == "pointinterval") {
    study.density.plot <- study.density.plot +
      ggdist::stat_pointinterval(
        ggplot2::aes(x = x_studies, y = Author),
        data = posterior.draws |> dplyr::filter(if (isTRUE(subgroup)) {Author_pooled != "Pooled Effect" & Author_pooled != "Overall Effect"}
                                                else {Author != "Pooled Effect"}),
        linewidth = 1,
        size = 1,
        color = color_pointinterval,
        limits = calc_xlim
      )
  }
  
  return(study.density.plot)
}

#' Internal function to summarise data for forest plot
#'
#' @noRd
forest.data.summary_fn <- function(spread_df,
                                   data,
                                   measure,
                                   sort_studies_by = sort_studies_by,
                                   subgroup = FALSE) {
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
  
  # Select only needed join columns
  join_vars <- unique(c("Author", "Year", props$data_cols))
  forest.data.summary <- forest.data |>
    dplyr::left_join(data |> dplyr::select(dplyr::any_of(join_vars), yi, vi, D1:Overall), by = "Author") |>
    dplyr::left_join(tau.summary, by = "Author", suffix = c("", "_sd")) |> 
    sort.studies.fn(sort_studies_by)
  
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
                                 paste0("τ = ", sprintf("%.2f", forest.data.summary$sd_Author__Intercept),
                                        " [", sprintf("%.2f", forest.data.summary$.lower_sd), ", ", sprintf("%.2f", forest.data.summary$.upper_sd), "]"),
                                 unweighted_effect))
  
  # Add study group summary columns depending on measure type
  if (measure %in% c("MD", "SMD")) {
    forest.data.summary <- forest.data.summary |>
      dplyr::mutate(
        N_int = dplyr::case_when(
          Author == "Pooled Effect" ~ as.character(sum(N_Intervention, na.rm = TRUE)),
          TRUE ~ as.character(N_Intervention)),
        int_mean_sd = dplyr::case_when(
          Author != "Pooled Effect" ~ paste0(sprintf("%.2f", Mean_Intervention), " (", sprintf("%.2f", SD_Intervention), ")"),
          TRUE ~ NA_character_),
        N_ctrl = dplyr::case_when(
          Author == "Pooled Effect" ~ as.character(sum(N_Control, na.rm = TRUE)),
          TRUE ~ as.character(N_Control)),
        ctrl_mean_sd = dplyr::case_when(
          Author != "Pooled Effect" ~ paste0(sprintf("%.2f", Mean_Control), " (", sprintf("%.2f", SD_Control), ")"),
          TRUE ~ NA_character_)
      )
  } else {
    forest.data.summary <- forest.data.summary |>
      dplyr::mutate(
        control_outcome_frac = dplyr::case_when(
          Author == "Pooled Effect" ~ paste0(
            sum(Event_Control[Author != "Pooled Effect"], na.rm = TRUE), "/",
            sum(N_Control[Author != "Pooled Effect"], na.rm = TRUE)),
          TRUE ~ paste0(Event_Control, "/", N_Control)),
        int_outcome_frac = dplyr::case_when(
          Author == "Pooled Effect" ~ paste0(
            sum(Event_Intervention[Author != "Pooled Effect"], na.rm = TRUE), "/",
            sum(N_Intervention[Author != "Pooled Effect"], na.rm = TRUE)),
          TRUE ~ paste0(Event_Intervention, "/", N_Intervention)))
  }
  
  # Remove rows where Author == "No Pooled Effect"
  forest.data.summary <- forest.data.summary |>
    dplyr::filter(Author != "No Pooled Effect")
  
  return(forest.data.summary)
}



#' Internal function to create table for left side of forest plot
#'
#' @noRd
forest.table.left_fn <- function(forest.data.summary,
                                 subgroup = FALSE,
                                 label_control = "Control",
                                 label_intervention = "Intervention",
                                 measure = "OR",
                                 font = NULL) {
  
  is_continuous <- measure %in% c("MD", "SMD")
  
  # Define column labels based on effect type
  if (is_continuous) {
    control_label <- paste(label_control, "\nMean (SD)")
    int_label <- paste(label_intervention, "\nMean (SD)")
  } else {
    control_label <- paste(label_control, "\n(Events/Total)")
    int_label <- paste(label_intervention, "\n(Events/Total)")
  }
  
  # Choose correct columns for data
  if (is_continuous) {
    selected_cols <- c("Author", "Year", "N_int", "int_mean_sd", "N_ctrl", "ctrl_mean_sd")
  } else {
    selected_cols <- c("Author", "Year", "int_outcome_frac", "control_outcome_frac")
  }
  if (isTRUE(subgroup)) {
    selected_cols <- c("Author", "Year", "Subgroup", selected_cols[!(selected_cols %in% c("Author", "Year"))])
  }
  
  df <- forest.data.summary |> dplyr::select(dplyr::any_of(selected_cols))
  
  df <- df |> dplyr::mutate(
    Author = dplyr::if_else(
      is.na(Year),
      Author,
      paste0(Author, " (", Year, ")")
    )
  )
  
  # Replace NA values with blanks for the 'Pooled Effect' row
  if (is_continuous) {
    df <- df |> dplyr::mutate(
      int_mean_sd = dplyr::if_else(Author == "Pooled Effect", NA_character_, int_mean_sd),
      ctrl_mean_sd = dplyr::if_else(Author == "Pooled Effect", NA_character_, ctrl_mean_sd)
    )
  }
  
  # Base gt table
  if (isFALSE(subgroup)) {
    forest.table.left <- df |>
      gt::gt() |>
      gt::cols_label(Author = "Study") |>
      gt::cols_align(align = "left") |>
      gt::tab_style(
        style = gt::cell_text(style = "italic", weight = "bold"),
        locations = gt::cells_body(rows = Author == "Pooled Effect")) 
  } else {
    forest.table.left <- df |>
      gt::gt(groupname_col = "Subgroup") |>
      gt::tab_options(row_group.font.weight = "bold") |>
      gt::cols_label(Author = gt::md("Subgroup/ \nStudy")) |>
      gt::cols_align(align = "left") |>
      gt::tab_style(
        gt::cell_text(color = "white"),
        locations = gt::cells_row_groups(groups = "Overall")) |> 
      gt::tab_style(
        style = gt::cell_text(weight = "bold", style = "italic"),
        locations = gt::cells_body(rows = Author == "Overall Effect")) |>
      gt::tab_style(
        style = gt::cell_text(style = "italic",  weight = "bold", color = "grey60"),
        locations = gt::cells_body(rows = Author == "Pooled Effect")) |>
      gt::tab_style(
        style = gt::cell_text(indent = gt::px(10)),
        locations = gt::cells_body(rows = Author != "Overall Effect",
                                   columns = Author))
  }
  
  # Column label mapping
  col_labels <- if (is_continuous) {
    c(N_int = "N", int_mean_sd = gt::md(int_label),
      N_ctrl = "N", ctrl_mean_sd = gt::md(control_label))
  } else {
    c(control_outcome_frac = gt::md(control_label),
      int_outcome_frac = gt::md(int_label))
  }
  
  # Final styling
  forest.table.left <- forest.table.left |>
    gt::cols_label(!!!col_labels) |>
    gt::tab_options(
      column_labels.font.weight = "bold",
      table.border.top.color = "white",
      table.border.bottom.color = "white",
      table.font.names = font) |>
    gt::cols_hide(columns = Year) |>
    gt::opt_table_lines(extent = "none") |> 
    gt::tab_style(
      style = gt::cell_fill(color = "grey95"),
      locations = gt::cells_body(rows = Author %in% c("Pooled Effect", "Overall Effect"))
    )
  
  if (measure %in% c("MD", "SMD")) {
    forest.table.left <- forest.table.left |>
      gt::sub_missing(columns = c(int_mean_sd, ctrl_mean_sd), missing_text = "") 
  }
  
  return(forest.table.left)
}


  #' Internal function to create table on right of forest plot (including RoB)
  #'
  #' @noRd  
forest.table.right_fn <- function(df,
                                  subgroup = FALSE,
                                  measure = "OR",
                                  add_rob = FALSE,
                                  font = NULL) {  
  if (isFALSE(subgroup)) {
    forest.table.right <- df |>
      dplyr::select(Author, weighted_effect, unweighted_effect, D1:Overall) |>
      gt::gt() |>
      gt::tab_style(
        style = list(
          gt::cell_text(
            style = "italic",
            weight = "bold")),
        locations = gt::cells_body(
          rows = Author == "Pooled Effect"))
    
  } else if (isTRUE(subgroup)) {
    forest.table.right <- df |>
      dplyr::select(Author, Subgroup, weighted_effect, unweighted_effect, D1:Overall) |>
      gt::gt(groupname_col = "Subgroup") |>
      gt::tab_style(
        gt::cell_text(color = "white"),
        locations = gt::cells_row_groups(groups = gt::everything())) |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold", style = "italic"),
        locations = gt::cells_body(rows = Author == "Overall Effect")) |>
      gt::tab_style(
        style = list(
          gt::cell_text(
            style = "italic", 
            weight = "bold",
            color = "grey60")),
        locations = gt::cells_body(
          rows = Author =="Pooled Effect"))
  }  
  
  forest.table.right <- forest.table.right |>
    gt::tab_options(column_labels.font.weight = "bold") |>
    gt::cols_align(align = "right") |>
    gt::cols_label(.list = create_rob_labels(df, measure)) |>
    gt::tab_options(
      table.border.top.color = "white",
      table.border.bottom.color = "white",
      table.font.names = font) |>
    gt::opt_table_lines(extent = "none") |>
    gt::tab_style(
      style = gt::cell_fill(color = "grey95"),
      locations = gt::cells_body(
        rows = Author %in% c("Pooled Effect", "Overall Effect"),
        columns = c(weighted_effect, unweighted_effect)))
  
  if (add_rob == FALSE) {
    forest.table.right <- forest.table.right |>
      gt::cols_hide(columns = c(Author, D1:Overall))
    
  } else if (add_rob == TRUE) {
    forest.table.right <- forest.table.right |>
      gt::opt_table_lines(extent = "none") |>
      gt::sub_missing(columns = c(D1:Overall), missing_text = "") |>
      gt::cols_align(
        align = "center",
        columns = c(D1:Overall)) |>
      # Reduce horizontal padding only for columns with "Risk", "of", "Bias" text
      gt::tab_style(
        style = "padding-left:0px; padding-right:0px;",
        locations = gt::cells_body(columns = get_rob_text_columns(df))) |>
      gt::tab_style(
        style = "padding-left:0px; padding-right:0px;",
        locations = gt::cells_column_labels(columns = get_rob_text_columns(df))) |> 
      gt::text_case_match(
        "High" ~ fontawesome::fa(name = "circle-plus", fill = "#e32400", height = "1.1em"),
        "Low" ~ fontawesome::fa(name = "circle-minus", fill = "#77bb41", height = "1.1em"),
        "Some concerns" ~ fontawesome::fa(name = "circle-question", fill = "#f5ec00", height = "1.1em")) |>
      gt::tab_style(
        style = "padding-top:2px; padding-bottom:2px;",
        locations = gt::cells_body(columns = c(D1:Overall))) |>
      gt::tab_style(
        style = "padding-right:2px;",
        locations = gt::cells_body(
          columns = unweighted_effect, 
          rows = Author == "Pooled Effect")) |>
      gt::cols_hide(columns = Author)
  }  

  return(forest.table.right)
}

#' Internal function to get columns that contain Risk/of/Bias text
#'
#' @noRd
get_rob_text_columns <- function(df) {
  
  # Get the column names that start with D or are "Overall"
  rob_cols <- names(df)[grepl("^D\\d+$|^Overall$", names(df))]
  
  # Split "Risk of Bias" across the middle columns
  rob_text <- c("Risk", "of", "Bias")
  n_rob_cols <- length(rob_cols)
  
  # Calculate which columns get the text
  if (n_rob_cols >= 3) {
    start_pos <- ceiling((n_rob_cols - 3) / 2) + 1
    text_positions <- start_pos:(start_pos + 2)
  } else {
    text_positions <- 1:min(n_rob_cols, 3)
  }
  
  # Return only the columns that have the Risk/of/Bias text
  return(rob_cols[text_positions])
}

#' Internal function to create RoB labels based on how many domains are present
#'
#' @noRd
create_rob_labels <- function(df, measure) {
  
  # Get the column names that start with D or are "Overall"
  rob_cols <- names(df)[grepl("^D\\d+$|^Overall$", names(df))]
  
  # Create base labels
  base_labels <- list(
    weighted_effect = gt::md(paste("Shrinkage", measure, "\n", "[95% CrI]")),
    unweighted_effect = gt::md(paste("Observed", measure, "\n", "[95% CrI]"))
  )
  
  # Split "Risk of Bias" across the middle columns
  rob_text <- c("Risk", "of", "Bias")
  n_rob_cols <- length(rob_cols)
  
  # Calculate which columns get the text
  if (n_rob_cols >= 3) {
    start_pos <- ceiling((n_rob_cols - 3) / 2) + 1
    text_positions <- start_pos:(start_pos + 2)
  } else {
    text_positions <- 1:min(n_rob_cols, 3)
    rob_text <- rob_text[1:length(text_positions)]
  }
  
  # Create labels for ROB columns using purrr
  rob_labels <- purrr::imap(rob_cols, ~ {
    col_name <- .x
    position <- .y
    
    if (position %in% text_positions) {
      text_idx <- which(text_positions == position)
      gt::md(paste0("**", rob_text[text_idx], "**\n", col_name))
    } else {
      gt::md(paste0("&nbsp;\n", col_name))
    }
  }) |>
    purrr::set_names(rob_cols)
  
  # Combine base labels with ROB labels
  return(c(base_labels, rob_labels))
}

#' Internal function to join sections of forest plot
#'
#' @noRd
patchwork_fn <- function(table.left,
                         study.density.plot,
                         table.right,
                         add_rob_legend = FALSE,
                         rob_tool = "rob2",
                         plot_width = NULL,
                         title = NULL,
                         subtitle = NULL,
                         title_align = "center",
                         font = NULL){
  if (is.null(plot_width)){
    plot_width <- 4
  }
  
  # Create the base plot layout
  if (isFALSE(add_rob_legend)) {
    base_plot <- patchwork::wrap_table(table.left, space = "fixed") + 
      study.density.plot + 
      patchwork::wrap_table(table.right, space = "fixed") +
      patchwork::plot_layout(widths = grid::unit(c(-1, plot_width, -1), 
                                                 c("null", "cm", "null")))
  } else {
    rob_legend <- rob_legend_fn(rob_tool, font)
    # Include rob_legend to the right of table.right
    base_plot <- patchwork::wrap_table(table.left, space = "fixed") + 
      study.density.plot + 
      patchwork::wrap_table(table.right, space = "fixed") +
      patchwork::wrap_table(rob_legend, space = "fixed") +
      patchwork::plot_layout(widths = grid::unit(c(-1, plot_width, -1, -1), 
                                                 c("null", "cm", "null", "null")))
  }
  
  # Add title and/or subtitle if provided
  if (!is.null(title) || !is.null(subtitle)) {
    # Set alignment value based on title_align parameter
    hjust_val <- switch(title_align,
                        "left" = 0,
                        "center" = 0.5,
                        "centre" = 0.5,
                        "right" = 1,
                        0.5) # default to center if invalid input
    
    # Create theme elements
    title_theme <- ggplot2::theme()
    
    if (!is.null(title)) {
      title_theme <- title_theme + 
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = 16, 
            face = "bold", 
            hjust = hjust_val,
            margin = ggplot2::margin(b = ifelse(is.null(subtitle), 10, 5)),
            family = font))
    }
    
    if (!is.null(subtitle)) {
      title_theme <- title_theme + 
        ggplot2::theme(
          plot.subtitle = ggplot2::element_text(
            size = 14, 
            hjust = hjust_val,
            margin = ggplot2::margin(b = 10),
            color = "gray30",
            family = font))
    }
    
    base_plot <- base_plot + 
      patchwork::plot_annotation(
        title = title,
        subtitle = subtitle,
        theme = title_theme
      )
  }
  
  return(base_plot)
}

#' Internal function to create risk of bias legend table
#'
#' @noRd
rob_legend_fn <- function(rob_tool = c("rob2", "robins_i", "quadas2", "robins_e", "nos"),
                          font = NULL) {
  
  # Get domains for the specified tool
  domains <- get_rob_domains(rob_tool, full_text = TRUE)
  
  # Remove overall domain
  domains <- domains[names(domains) != "Overall"]
  
  # Create the data for the table
  legend_data <- data.frame(
    Code = paste0("(", names(domains), ")"),
    Description = unname(domains),
    stringsAsFactors = FALSE)
  
  # Add empty row for spacing
  legend_data <- rbind(
    legend_data,
    data.frame(Code = "", Description = "", stringsAsFactors = FALSE))
  
  # Add risk level rows with placeholder text
  risk_rows <- data.frame(
    Code = c(" ", "Low", "Some concerns", "High"),
    Description = c(" ", "Low risk", "Some concerns", "High risk"),
    stringsAsFactors = FALSE)
  
  legend_data <- rbind(legend_data, risk_rows)
  
  # Create the gt table
  rob_legend_table <- legend_data |>
    gt::gt() |>
    # Remove column headers
    gt::cols_label(
      Code = "",
      Description = "") |>
    # Add title
    gt::tab_header(
      title = gt::html("<b>Risk of bias legend<b>")) |>
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold", size = gt::px(16))),
      locations = gt::cells_title(groups = "title")) |> 
    # Style the domain codes (first n rows)
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold", size = gt::px(12))),
      locations = gt::cells_body(
        columns = Code,
        rows = 1:length(domains))) |>
    # Transform the risk level text to icons
    gt::text_transform(
      locations = gt::cells_body(
        columns = Code,
        rows = (length(domains) + 2):nrow(legend_data)),
      fn = function(x) {
        dplyr::case_when(
          x == "High" ~ as.character(fontawesome::fa(name = "circle-plus", fill = "#e32400", height = "1.5em")),
          x == "Low" ~ as.character(fontawesome::fa(name = "circle-minus", fill = "#77bb41", height = "1.5em")),
          x == "Some concerns" ~ as.character(fontawesome::fa(name = "circle-question", fill = "#f5ec00", height = "1.5em")),
          TRUE ~ x)
      }) |>
    # Style the risk level descriptions to be bold and italic
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold", style = "italic")),
      locations = gt::cells_body(
        columns = Description,
        rows = (length(domains) + 2):nrow(legend_data))) |>
    # Reduce space in between RoB and legend
    gt::tab_style(
      style = "padding-left:2px;",
      locations = gt::cells_body(
        columns = Code, 
        rows = everything())) |> 
    # Set column widths
    gt::cols_width(
      Code ~ gt::px(40),
      Description ~ gt::px(300)) |>
    # Table options
    gt::tab_options(
      table.border.top.color = "white",
      table.border.bottom.color = "white",
      heading.align = "left",
      heading.padding = gt::px(6),
      table.font.size = gt::px(12),
      table.font.names = font) |>
    gt::opt_table_lines(extent = "none")
  
  return(rob_legend_table)
}


#' Internal function to get RoB domains for RoB legend
#'
#' @noRd
get_rob_domains <- function(rob_tool = c("rob2", "robins_i", "quadas2", "robins_e"),
                            full_text = TRUE) {
  rob_tool <- match.arg(rob_tool)

  domains <- switch(rob_tool,
                    rob2 = c(
                      D1 = "Bias arising from the randomization process",
                      D2 = "Bias due to deviations from intended interventions",
                      D3 = "Bias due to missing outcome data",
                      D4 = "Bias in measurement of the outcome",
                      D5 = "Bias in selection of the reported result",
                      Overall = "Overall risk of bias"
                    ),
                    robins_i = c(
                      D1 = "Bias due to confounding",
                      D2 = "Bias in selection of participants into the study",
                      D3 = "Bias in classification of interventions",
                      D4 = "Bias due to deviations from intended interventions",
                      D5 = "Bias due to missing data",
                      D6 = "Bias in measurement of outcomes",
                      D7 = "Bias in selection of the reported result",
                      Overall = "Overall risk of bias"
                    ),
                    quadas2 = c(
                      D1 = "Patient selection",
                      D2 = "Index test(s)",
                      D3 = "Reference standard",
                      D4 = "Flow and timing",
                      Overall = "Overall risk of bias"
                    ),
                    robins_e = c(
                      D1 = "Bias due to confounding",
                      D2 = "Bias in selection of participants into the study",
                      D3 = "Bias in classification of exposures",
                      D4 = "Bias due to departures from intended exposures",
                      D5 = "Bias due to missing data",
                      D6 = "Bias in measurement of outcomes",
                      D7 = "Bias in selection of the reported result",
                      Overall = "Overall risk of bias"
                    )
  )

  if (isTRUE(full_text)){ 
    return(domains)
    } else {
      return(names(domains))
    }
}
