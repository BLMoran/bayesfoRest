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

forest.data_fn <- function(data,
                           model,
                           subgroup = FALSE,
                           sort_studies_by = "author",
                           custom_group_order = NULL) {

  if (subgroup == FALSE && "Subgroup" %in% names(data)) {
    data <- data |> dplyr::select(-Subgroup)
  }

  if (subgroup == FALSE) {
    # No subgroups, evaluate as a whole
    study.draws <- tidybayes::spread_draws(model, r_Author[Author, ], b_Intercept) |>
      dplyr::mutate(b_Intercept = r_Author + b_Intercept)

    pooled.draws <- tidybayes::spread_draws(model, b_Intercept, sd_Author__Intercept) |>
      dplyr::mutate(Author = "Pooled Effect")

    effect.draws <- dplyr::bind_rows(study.draws, pooled.draws) |>
      dplyr::ungroup() |>
      dplyr::left_join(dplyr::select(data, Author, Year, yi, vi), by = "Author") |>
      sort.studies.fn(sort_studies_by)

  } else {
    # With subgroups
    subgroup_df <- data |>
      dplyr::group_by(Subgroup) |>
      tidyr::nest() |>
      dplyr::mutate(
        subgroup_model = purrr::map(data, ~ update(model, newdata = .x))
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
        effect_draws = purrr::map2(subgroup_model, data, ~ {
          study <- tidybayes::spread_draws(.x, r_Author[Author, ], b_Intercept) |>
            dplyr::mutate(b_Intercept = r_Author + b_Intercept)

          pooled <- tidybayes::spread_draws(.x, b_Intercept, sd_Author__Intercept) |>
            dplyr::mutate(Author = "Pooled Effect")

          combined <- dplyr::bind_rows(study, pooled) |>
            dplyr::ungroup() |>
            dplyr::left_join(dplyr::select(.y, Author, Year, yi, vi), by = "Author")|>
            sort.studies.fn(sort_studies_by="year")

          return(combined)
        })
      ) |>
      tidyr::unnest(effect_draws) |>
      dplyr::select(-data, -subgroup_model)

    effect.draws <- dplyr::bind_rows(study.effect.draws, overall.effect.draws)

    # Custom group order for Subgroup column
    if (!is.null(custom_group_order)) {
      effect.draws <- effect.draws |>
        dplyr::mutate(Subgroup = factor(Subgroup, levels = custom_group_order)) |>
        dplyr::arrange(Subgroup) |>
        dplyr::mutate(Subgroup =  dplyr::case_when(
          is.na(Subgroup) & Author == "Overall Effect" ~ "Overall",
          .default = Subgroup
        ))
    }
  }

  return(effect.draws)
}


study.density.plot_fn  <- function(df, effectsize = "OR",
                                   subgroups = FALSE,
                                   custom_group_order = NULL,
                                   color_study_posterior = "dodgerblue",
                                   color_study_posterior_outline = "blue",
                                   color_overall_posterior = "blue",
                                   color_shrinkage_outline = "purple",
                                   color_pointinterval = "purple",
                                   color_shrinkage_fill = NA,
                                   label_control = "Control",
                                   label_intervention = "Intervention",
                                   shrinkage_output = "density"){


  # Optimise data structure: separate study-level data from posterior draws
  if (isTRUE(subgroups)){
    study.effects <- df |> dplyr::distinct(Author_ordered, yi, vi) |>
      dplyr::mutate(
        Author_ordered = factor(Author_ordered, levels = as.character(1:max(Author_ordered, na.rm = TRUE))),
        Author = fct_rev(Author_ordered))

    posterior.draws <- df |>
      dplyr::mutate(
        Author_pooled = Author,
        Author_ordered = factor(Author_ordered, levels = as.character(1:max(Author_ordered, na.rm = TRUE))),
        Author = fct_rev(Author_ordered))


  } else if (isFALSE(subgroups)){
    study.effects <- df |> dplyr::distinct(Author, yi, vi)
    posterior.draws <- df
  }

  # Calculate x.min based on effectsize
  x.min <- if (effectsize == "OR") {
    0.1
  } else {
    floor(min(exp(df$yi - 1.96 * sqrt(df$vi)), na.rm = TRUE))
  }

  # Calculate x.max
  x.max <- ceiling(max(exp(df$yi + 1.96 * sqrt(df$vi)), na.rm = TRUE))

  # Set xlim to given value or calculated range
  calc_xlim <- c(x.min, x.max)

  study.density.plot <-
    ggplot(aes(y = Author), data = posterior.draws) +
    ggdist::stat_slab(aes(xdist = exp(distributional::dist_normal(mean = yi, sd = sqrt(vi)))),
                      slab_linewidth = 0.5, alpha = 0.7, limits = calc_xlim,
                      data = study.effects, colour = color_study_posterior_outline, fill = color_study_posterior) +
    ggdist::stat_slab(aes(x = exp(b_Intercept), y = Author),
                      data = posterior.draws |> dplyr::filter(if (isTRUE(subgroups)) {Author_pooled == "Pooled Effect"|Author_pooled == "Overall Effect"}
                                                           else {Author == "Pooled Effect"}),
                      fill = color_overall_posterior, scale = 0.8) +
    guides(alpha = "none", fill = "none")+
    # Add vertical lines effect and CI
    geom_vline(xintercept = exp(fixef(model)[1, 1]), color = "grey60", linewidth = 1) +
    geom_vline(xintercept = exp(fixef(model)[1, 3:4]), color = "grey60", linetype = 2) +
    geom_vline(xintercept = 1, color = "black", linewidth = 1) +
    scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 7), expand = c(0, 0), limits = calc_xlim)+
    coord_cartesian(xlim= calc_xlim, clip = "off") +
    theme_light() +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          axis.title.x = element_text(vjust=-0.5), panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = margin(0,0,0,0), panel.border = element_blank(),
          axis.line.x.top = element_line(color = "grey60", linewidth = 0.75),
          axis.line.x.bottom = element_line(color = "black", linewidth = 0.75),
          axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank(),
          axis.text.x.bottom = element_text(colour = "black")) +
    guides(x.sec = "axis", y.sec = "axis") +
    annotation_custom(grid::textGrob(
      label = paste(" Favours\n", label_intervention),
      x = unit(-0.01, "npc"), y = unit(1.02, "npc") , just = c("left", "bottom"),
      gp = grid::gpar( col = "grey30")),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
    annotation_custom(grid::textGrob(
      label = paste(" Favours\n", label_control),
      x = unit(1, "npc"), y = unit(1.02, "npc") , just = c("right", "bottom"),
      gp = grid::gpar(col = "grey30")),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
    labs(x="Odds ratio (log scale)") +
    ylab(NULL) +
    geom_hline(yintercept = 1, color = "black", linewidth = 0.75)

  if(isFALSE(subgroups)){
    study.density.plot <- study.density.plot + scale_y_discrete(expand = c(0, 0), limits = rev)
  } else if (isTRUE(subgroups)){
    study.density.plot <- study.density.plot + scale_y_discrete(expand = c(0, 0))
  }

  if(is.null(shrinkage_output))"density" else shrinkage_output

  if (shrinkage_output == "density") {
    study.density.plot <- study.density.plot +
      ggdist::stat_slab(
        aes(x = exp(b_Intercept), y = Author),
        data = posterior.draws |> dplyr::filter(if (isTRUE(subgroups)) {Author_pooled != "Pooled Effect" & Author_pooled != "Overall Effect"}
                                                else {Author != "Pooled Effect"}),
        linewidth = 0.5,
        scale = 0.8,
        color = color_shrinkage_outline,
        fill = color_shrinkage_fill,
        limits = calc_xlim
      )

  } else if (shrinkage_output == "pointinterval") {
    study.density.plot <- study.density.plot +
      ggdist::stat_pointinterval(
        aes(x = exp(b_Intercept), y = Author),
        data = posterior.draws |> dplyr::filter(if (isTRUE(subgroups)) {Author_pooled != "Pooled Effect" & Author_pooled != "Overall Effect"}
                                                else {Author != "Pooled Effect"}),
        linewidth = 1,
        size = 1,
        color = color_pointinterval,
        limits = calc_xlim
      )
  }

  return(study.density.plot)
}


forest.data.summary_fn <- function(spread_df,
                                   data,
                                   sort_studies_by = "author",
                                   subgroups = FALSE) {

  # Study Summaries
  forest.data <- spread_df |>
    dplyr::group_by(Author) |>
    tidybayes::median_hdi(b_Intercept)

  # Tau Summary
  tau.summary <- spread_df |>
    dplyr::group_by(Author) |>
    tidybayes::median_hdi(sd_Author__Intercept)

  forest.data.summary <- forest.data |>
    dplyr::left_join(data |> dplyr::select(Author, Year, yi, vi, N_Total, N_Control, N_Intervention, Outcome_Control_Yes, Outcome_Intervention_Yes, D1:Overall), by = "Author") |>
    dplyr::mutate(
      control_outcome_frac = ifelse(Author == "Pooled Effect",
                                    paste0(sum(Outcome_Control_Yes[Author != "Pooled Effect"], na.rm = TRUE),
                                           "/", sum(N_Control[Author != "Pooled Effect"], na.rm = TRUE)),
                                    paste0(Outcome_Control_Yes, "/", N_Control)),
      int_outcome_frac = ifelse(Author == "Pooled Effect",
                                paste0(sum(Outcome_Intervention_Yes[Author != "Pooled Effect"], na.rm = TRUE),
                                       "/", sum(N_Intervention[Author != "Pooled Effect"], na.rm = TRUE)),
                                paste0(Outcome_Intervention_Yes, "/", N_Intervention)),
      weighted_effect = paste0(
        sprintf("%.2f", exp(b_Intercept)), " [",
        sprintf("%.2f", exp(.lower)), ", ",
        sprintf("%.2f", exp(.upper)), "]"
      ),
      unweighted_effect = paste0(
        sprintf("%.2f", exp(yi)), " [",
        sprintf("%.2f", exp(yi - 1.96 * sqrt(vi))), ", ",
        sprintf("%.2f", exp(yi + 1.96 * sqrt(vi))), "]"
      )
    ) |>
    dplyr::mutate(unweighted_effect = ifelse(unweighted_effect == "NA [NA, NA]",
                                             paste0("τ = ", sprintf('%.2f', tau.summary$sd_Author__Intercept),
                                                    " [", sprintf('%.2f', tau.summary$.lower), ", ",
                                                    sprintf('%.2f', tau.summary$.upper), "]"),
                                             unweighted_effect)) |>
    sort.studies.fn(sort_studies_by)

  if (isFALSE(subgroups) && "Subgroup" %in% names(forest.data.summary)) {
    forest.data.summary <- forest.data.summary |>
      dplyr::select(-Subgroup)
  }

  return(forest.data.summary)
}


# Function to create table to left of density plots
forest.table.left_fn <- function(forest.data.summary,
                                 subgroups = FALSE,
                                 label_control = "Control",
                                 label_intervention = "Intervention",
                                 effectsize = "OR") {

  control_label <- if (effectsize %in% c("OR", "HR")) {
    paste(label_control, "\n", "(Events/Total)")
  } else {
    paste(label_control, "\n", "Mean (SD)")
  }

  int_label <- if (effectsize %in% c("OR", "HR")) {
    paste(label_intervention, "\n", "(Events/Total)")
  } else {
    paste(label_intervention, "\n", "Mean (SD)")
  }

  if(isFALSE(subgroups)){
    subgroup.forest.table.left <- forest.data.summary |>
      dplyr::select(Author, Year, int_outcome_frac, control_outcome_frac) |>
      dplyr::mutate(
        Author = dplyr::if_else(
          is.na(Year),
          Author,
          paste0(Author, " (", Year, ")")
        )
      ) |>
      gt::gt() |>
      gt::cols_label(Author = "Study")

  } else if (isTRUE(subgroups)){
    subgroup.forest.table.left <- forest.data.summary |>
      dplyr::select(Author, Year, Subgroup, int_outcome_frac, control_outcome_frac) |>
      dplyr::mutate(
        Author = dplyr::if_else(
          is.na(Year),
          Author,
          paste0(Author, " (", Year, ")")
        )
      ) |>
      gt::gt(groupname_col = "Subgroup") |>
      gt::tab_options(row_group.font.weight = "bold") |>
      gt::cols_label(Author = gt::md("Subgroup/ \nStudy")) |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold", style = "italic"),
        locations = gt::cells_body(rows = Author == "Overall Effect")
      )
  }

  subgroup.forest.table.left <- subgroup.forest.table.left |>
    gt::tab_options(column_labels.font.weight = "bold") |>
    gt::cols_align(align = "left") |>
    gt::cols_label(
      control_outcome_frac = gt::md(control_label),
      int_outcome_frac = gt::md(int_label)
    ) |>
    gt::tab_style(
      style = gt::cell_text(style = "italic", color = "grey60"),
      locations = gt::cells_body(rows = Author == "Pooled Effect")
    ) |>
    gt::tab_options(
      table.border.top.color = "white",
      table.border.bottom.color = "white"
    ) |>
    gt::cols_hide(columns = Year) |>
    gt::opt_table_lines(extent = "none")
}


forest.table.right_fn <- function(df,
                                  subgroups = FALSE,
                                  effectsize = "OR",
                                  add_rob=FALSE,
                                  rob_legend = FALSE,
                                  rob_tool = "rob2"){

  if(isFALSE(subgroups)){
    forest.table.right <- df |>
      dplyr::select(Author, weighted_effect, unweighted_effect, D1:Overall) |>
      gt::gt()

  } else if (isTRUE(subgroups)){
    forest.table.right <- df |>
      dplyr::select(Author, Subgroup, weighted_effect, unweighted_effect, D1:Overall) |>
      gt::gt(groupname_col = "Subgroup") |>
      gt::tab_style(
        gt::cell_text(color = "white"),
        locations = gt::cells_row_groups(groups = everything())) |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold", style = "italic"),
        locations = gt::cells_body(rows = Author == "Overall Effect")
      )
  }

    forest.table.right <- forest.table.right |>
      gt::tab_options(column_labels.font.weight = "bold") |>
      gt::cols_align(align = "left") |>
      gt::cols_label(weighted_effect =gt:: md(paste("Shrinkage", effectsize, "\n", "[95% CrI]")),
                 unweighted_effect = gt::md(paste("Observed", effectsize, "\n", "[95% CrI]"))) |>
      gt::tab_style(
        style = list(
          gt::cell_text(style = "italic",
                    color = "grey60")),
        locations = gt::cells_body(
          rows = Author == "Pooled Effect")) |>
      gt::tab_options(table.border.top.color = "white",
                  table.border.bottom.color = "white") |>
      gt::opt_table_lines(extent = "none")

  if (add_rob == FALSE){
    forest.table.right <- forest.table.right |>
      gt::cols_hide(columns = c(Author, D1:Overall))

  } else if (add_rob == TRUE) {
    forest.table.right <-  forest.table.right |>
      gt::opt_table_lines(extent = "none") |>
      gt::sub_missing(columns = c(D1:Overall), missing_text = "") |>
      gt::cols_align(align = "center",
                 columns = c(D1:Overall)) |>
      gt::tab_spanner(
        label = gt::md("**Risk of Bias**"),
        columns = c(D1:Overall),
        id = "rob") |>
      gt::tab_style(
        style = "padding-top:2px;",
        locations = gt::cells_column_labels(columns = c(D1:Overall))) |>
      gt::tab_style(
        style = list(
          "padding-bottom:2px;",
          gt::cell_text(weight = "bold")),
        locations = gt::cells_column_spanners(spanners = "rob")) |>
      gt::text_case_match(
        "High" ~ fontawesome::fa(name="circle-plus", fill = "#e32400", height = "1.1em"),
        "Low" ~ fontawesome::fa(name="circle-minus", fill= "#77bb41", height = "1.1em"),
        "Unsure" ~ fontawesome::fa(name="circle-question", fill="#f5ec00", height = "1.1em")) |>
      gt::tab_style(
        style = "padding-top:2px; padding-bottom:2px;",
        locations = gt::cells_body(columns = c(D1:Overall))) |>
      gt::tab_style(
        style = "padding-right:2px;",
        locations = gt::cells_body(columns = unweighted_effect, rows = Author == "Pooled Effect")) |>
      gt::cols_hide(columns = Author)
  }

    if (isTRUE(add_rob) && isTRUE(rob_legend)) {
      forest.table.right <- forest.table.right |>
        add_rob_legend(rob_tool = rob_tool)
    }

  return(forest.table.right)
}



patchwork_fn <- function(table.left,
                         study.density.plot,
                         table.right,
                         plot_width = NULL){

  if (is.null(plot_width)){
    plot_width <- 4
  } else{
    plot_width <- plot_width
  }

  patchwork::wrap_table(table.left, space = "fixed") + study.density.plot + patchwork::wrap_table(table.right, space = "fixed") +
    patchwork::plot_layout(widths = unit(c(-1, plot_width,-1), c("null","cm","null")))
}

add_rob_legend <- function(.gt_tbl,
                           rob_tool = "rob2") {

  rob_domains <- get_rob_domains(rob_tool, full_text = TRUE)
  rob_labels <- names(rob_domains)

  rob_legend_text <- paste0(
    purrr::map2_chr(rob_labels, rob_domains, ~ paste0(.x, ": ", .y)),
    collapse = "<br>")

  .gt_tbl |>
    gt::tab_source_note(
      gt::md(paste0("**Risk of Bias Domains**<br>", rob_legend_text))
    )
}


get_rob_domains <- function(rob_tool = c("rob2", "robins_i", "quadas2", "robins_e", "nos"),
                            full_text = TRUE) {
  rob_tool <- match.arg(tolower(rob_tool))  # FIXED this line

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
                    ),
                    nos = c(
                      D1 = "Selection of study groups",
                      D2 = "Comparability of groups",
                      D3 = "Outcome (cohort) or exposure (case-control)",
                      Overall = "Overall study quality (via stars)"
                    )
  )

  if (full_text) return(domains) else return(names(domains))
}

