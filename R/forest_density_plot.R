#' Internal function to generate plot section of forest plot
#'
#' @description
#' Creates the density plot component of the forest plot. Uses the canonical
#' draws format (theta_study, eta, tau) and a pre-computed pooled_effect
#' summary instead of calling brms directly.
#'
#' @param df Draws data frame with canonical columns (theta_study, etc.).
#' @param pooled_effect Named numeric vector with estimate, lower, upper
#'   (on model scale, from extract_pooled_summary()).
#' @param measure Effect measure string.
#' @param effect_type "random_effect" or "common_effect".
#' @param subgroup Logical.
#' @param subgroup_order Subgroup ordering vector.
#' @param plot_output "density" or "boxplot". Boxplots are horizontal: the box
#'   spans the interquartile range, the line marks the median, and the whiskers
#'   span the 95% interval.
#' @param shrinkage_output "density" or "pointinterval".
#' @param ... Additional colour and layout arguments.
#'
#' @noRd
study.density.plot_fn  <- function(df,
                                   pooled_effect,
                                   measure = measure,
                                   effect_type = "random_effect",
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
                                   split_color_by_null = FALSE,
                                   color_favours_control = "firebrick",
                                   color_favours_intervention = "dodgerblue",
                                   label_control = "Control",
                                   label_intervention = "Intervention",
                                   add_arm_labels = TRUE,
                                   reverse_arm_labels = FALSE,
                                   plot_output = "density",
                                   shrinkage_output = "density",
                                   add_pred = FALSE,
                                   pred_output = "density",
                                   color_pred_posterior = "orange",
                                   color_pred_outline = "darkorange2",
                                   color_pred_pointinterval = "orange",
                                   xlim = NULL,
                                   x_breaks = NULL,
                                   x_label = NULL,
                                   null_value = NULL,
                                   add_rope = FALSE,
                                   rope_value = NULL,
                                   rope_color = "grey50",
                                   font = NULL){

  # Filter out "No Pooled Effect" rows at the beginning
  df <- df |> dplyr::filter(Author != "No Pooled Effect")

  # Apply color palette if provided
  if (!is.null(color_palette)) {
    palette_colors <- paletteer::paletteer_d(color_palette, n = 6)
    color_study_posterior <- palette_colors[2]
    color_study_posterior_outline <- palette_colors[1]
    color_pooled_posterior <- palette_colors[3]
    color_overall_posterior <- palette_colors[5]
    color_shrinkage_outline <- palette_colors[6]
    color_pointinterval <- palette_colors[6]
    if (!is.na(color_shrinkage_fill)) {
      color_shrinkage_fill <- scales::alpha(palette_colors[6], 0.3)
    }
  }

  # Get effect size properties
  props <- get_measure_properties(measure)

  # Determine x_breaks to use: custom if provided, otherwise use measure defaults
  breaks <- if (!is.null(x_breaks)) x_breaks else ggplot2::waiver()

  # Determine null line to use: custom if provided, otherwise use measure default
  null_value <- if (!is.null(null_value)) null_value else props$null_value

  # Process ROPE values if provided
  if (isTRUE(add_rope)) {
    if (is.null(rope_value)) {
      if (measure %in% c("OR", "HR", "RR", "IRR")) {
        rope_range <- c(0.9, 1.1)
      } else if (measure == "SMD") {
        rope_range <- c(-0.1, 0.1)
      } else if (measure == "MD") {
        stop("For mean differences, rope_value must be specified as it depends on the measurement scale")
      }
    } else {
      if (length(rope_value) == 1) {
        rope_range <- c(null_value - rope_value, null_value + rope_value)
      } else if (length(rope_value) == 2) {
        rope_range <- sort(rope_value)
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

  # Filter to study rows with valid yi/vi for axis limit calculation
  study_yi <- df$yi[!is.na(df$yi) & !is.na(df$vi)]
  study_vi <- df$vi[!is.na(df$yi) & !is.na(df$vi)]

  if (isTRUE(props$log_scale)){
    x.min <- 0.1
    x.max <- if (length(study_yi) > 0) ceiling(max(exp(study_yi + 1.96 * sqrt(study_vi)))) else 10
  } else if (isFALSE(props$log_scale)){
    x.min <- if (length(study_yi) > 0) floor(min(study_yi - 1.96 * sqrt(study_vi))) else -2
    x.max <- if (length(study_yi) > 0) ceiling(max(study_yi + 1.96 * sqrt(study_vi))) else 2
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
        x_studies = exp(theta_study))

    study.effects <- study.effects |>
      dplyr::mutate(
        xdist = exp(distributional::dist_normal(mean = yi, sd = sqrt(vi))))
  } else {
    posterior.draws <- posterior.draws |>
      dplyr::mutate(
        x_studies = theta_study)

    study.effects <- study.effects |>
      dplyr::mutate(
        xdist = distributional::dist_normal(mean = yi, sd = sqrt(vi)))
  }

  box_layers <- if (plot_output == "boxplot") {
    pooled_rows <- if (isTRUE(subgroup)) {
      posterior.draws |> dplyr::filter(Author_pooled %in% c("Pooled Effect", "Overall Effect"))
    } else {
      posterior.draws |> dplyr::filter(Author == "Pooled Effect")
    }
    box_df <- dplyr::bind_rows(
      study.effects |>
        dplyr::filter(!is.na(yi), !is.na(vi)) |>
        box_stats_normal() |>
        dplyr::mutate(.fill = color_study_posterior),
      pooled_rows |>
        box_stats_draws(overall = if (isTRUE(subgroup)) Author_pooled == "Overall Effect" else FALSE) |>
        dplyr::mutate(.fill = dplyr::if_else(.overall, color_overall_posterior, color_pooled_posterior)) |>
        dplyr::select(-.overall)) |>
      box_stats_scale(log_scale = isTRUE(props$log_scale), xlim = calc_xlim) |>
      dplyr::mutate(.above_null = xmiddle > null_value)

    if (isTRUE(split_color_by_null)) {
      ggplot2::geom_boxplot(
        ggplot2::aes(y = Author, xmin = xmin, xlower = xlower, xmiddle = xmiddle,
                     xupper = xupper, xmax = xmax, fill = .above_null),
        data = box_df, stat = "identity", orientation = "y", width = 0.35,
        position = ggplot2::position_nudge(y = 0.3), alpha = 0.7,
        colour = color_study_posterior_outline, linewidth = 0.5, inherit.aes = FALSE)
    } else {
      list(
        ggplot2::geom_boxplot(
          ggplot2::aes(y = Author, xmin = xmin, xlower = xlower, xmiddle = xmiddle,
                       xupper = xupper, xmax = xmax, fill = .fill),
          data = box_df, stat = "identity", orientation = "y", width = 0.35,
          position = ggplot2::position_nudge(y = 0.3), alpha = 0.7,
          colour = color_study_posterior_outline, linewidth = 0.5, inherit.aes = FALSE),
        ggplot2::scale_fill_identity())
    }
  }

  # Density slabs extend each row's range to n + 0.9; match it for boxplots
  # so rows line up with the tables
  box_ylim <- if (plot_output == "boxplot") {
    c(1, dplyr::n_distinct(posterior.draws$Author) + 0.9)
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
    {if (plot_output != "density") NULL else if (isTRUE(split_color_by_null)) {
      list(
        ggdist::stat_slab(
          ggplot2::aes(xdist = xdist, fill = ggplot2::after_stat(x > null_value)),
          slab_linewidth = 0.5, alpha = 0.7, limits = calc_xlim, height = 0.9,
          data = study.effects, colour = color_study_posterior_outline))
    } else {
      ggdist::stat_slab(
        ggplot2::aes(xdist = xdist),
        slab_linewidth = 0.5, alpha = 0.7, limits = calc_xlim, height = 0.9,
        data = study.effects, colour = color_study_posterior_outline,
        fill = color_study_posterior)
    }} +
    {if (plot_output != "density") NULL else if (isTRUE(split_color_by_null)) {
      list(
        ggdist::stat_slab(
          ggplot2::aes(
            x = x_studies, y = Author, fill = ggplot2::after_stat(x > null_value)),
          data = posterior.draws |>
            dplyr::filter(if (isTRUE(subgroup)) {
              Author_pooled == "Pooled Effect"
            } else {
              Author == "Pooled Effect"
            }),
          height = 0.9, normalize = "groups", colour = color_study_posterior_outline,
          limits = calc_xlim))
    } else {
      ggdist::stat_slab(
        ggplot2::aes(x = x_studies, y = Author),
        data = posterior.draws |>
          dplyr::filter(if (isTRUE(subgroup)) {
            Author_pooled == "Pooled Effect"
          } else {
            Author == "Pooled Effect"
          }),
        fill = color_pooled_posterior, height = 0.9, normalize = "groups",
        limits = calc_xlim)
    }} +
    # Add overall effect slab when subgroup = TRUE
    {if (isTRUE(subgroup) && plot_output == "density") {
      if (isTRUE(split_color_by_null)) {
        list(
          ggdist::stat_slab(
            ggplot2::aes(
              x = x_studies, y = Author, fill = ggplot2::after_stat(x > null_value)),
            data = posterior.draws |> dplyr::filter(Author_pooled == "Overall Effect"),
            height = 0.9, normalize = "groups", colour = color_study_posterior_outline))
      } else {
        ggdist::stat_slab(
          ggplot2::aes(x = x_studies, y = Author),
          data = posterior.draws |> dplyr::filter(Author_pooled == "Overall Effect"),
          fill = color_overall_posterior, height = 0.9, normalize = "groups")
      }
    }} +
    {if (plot_output == "boxplot") box_layers} +
    {if (isTRUE(split_color_by_null)) {
      ggplot2::scale_fill_manual(
        values = c(
          "FALSE" = color_favours_intervention,
          "TRUE"  = color_favours_control),
        guide = "none")
    }} +
    # Add null value
    ggplot2::geom_vline(xintercept = null_value, color = "black", linewidth = 1) +
    ggplot2::coord_cartesian(xlim = calc_xlim, ylim = box_ylim, clip = "off") +
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
    {if (isTRUE(add_arm_labels)) {
      # Swap labels if reversed
      right_label <- if (isTRUE(reverse_arm_labels)) label_intervention else label_control
      left_label  <- if (isTRUE(reverse_arm_labels)) label_control else label_intervention
      list(
        ggplot2::annotation_custom(grid::textGrob(
          label = paste(" Favours\n", right_label),
          x = grid::unit(1, "npc"), y = grid::unit(1.02, "npc"), just = c("right", "bottom"),
          gp = grid::gpar(col = "grey30", fontsize = 10, fontfamily = font)),
          xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
        ggplot2::annotation_custom(grid::textGrob(
          label = paste(" Favours\n", left_label),
          x = grid::unit(0, "npc"), y = grid::unit(1.02, "npc"), just = c("left", "bottom"),
          gp = grid::gpar(col = "grey30", fontsize = 10, fontfamily = font)),
          xmin = calc_xlim[1] - 0.01, xmax = calc_xlim[2], ymin = -Inf, ymax = Inf)
      )
    }} +
    ggplot2::labs(x = x_label %||% props$x_label) +
    ggplot2::ylab(NULL) +
    ggplot2::geom_hline(yintercept = 1, color = "black", linewidth = 0.75)

  # Pooled effect reference lines
  if (isTRUE(props$log_scale)){
    study.density.plot <- study.density.plot +
      ggplot2::scale_x_log10(breaks = breaks, expand = c(0, 0), limits = calc_xlim) +
      ggplot2::geom_vline(xintercept = exp(pooled_effect["estimate"]), color = "grey60", linewidth = 1) +
      ggplot2::geom_vline(xintercept = exp(pooled_effect[c("lower", "upper")]), color = "grey60", linetype = 2)
  } else {
    study.density.plot <- study.density.plot +
      ggplot2::scale_x_continuous(breaks = breaks, expand = c(0, 0), limits = calc_xlim) +
      ggplot2::geom_vline(xintercept = pooled_effect["estimate"], color = "grey60", linewidth = 1) +
      ggplot2::geom_vline(xintercept = pooled_effect[c("lower", "upper")], color = "grey60", linetype = 2)
  }

  if(isFALSE(subgroup)){
    study.density.plot <- study.density.plot + ggplot2::scale_y_discrete(expand = c(0, 0), limits = rev)
  } else if (isTRUE(subgroup)){
    # Fix the rows explicitly so empty subgroup header rows are kept
    study.density.plot <- study.density.plot +
      ggplot2::scale_y_discrete(expand = c(0, 0), limits = levels(posterior.draws$Author))
  }

  if(is.null(shrinkage_output))"density" else shrinkage_output

  # Shrinkage layer — only for random-effects models
  if (effect_type == "random_effect") {
    if (shrinkage_output == "density") {
      study.density.plot <- study.density.plot +
        ggdist::stat_slab(
          ggplot2::aes(x = x_studies, y = Author),
          data = posterior.draws |> dplyr::filter(if (isTRUE(subgroup)) {Author_pooled != "Pooled Effect" & Author_pooled != "Overall Effect" & Author_pooled != "Prediction"}
                                                  else {Author != "Pooled Effect" & Author != "Prediction"}),
          linewidth = 0.5, scale = 0.6, normalize = "groups",
          color = color_shrinkage_outline, fill = color_shrinkage_fill, limits = calc_xlim)
    } else if (shrinkage_output == "pointinterval") {
      study.density.plot <- study.density.plot +
        ggdist::stat_pointinterval(
          ggplot2::aes(x = x_studies, y = Author),
          data = posterior.draws |> dplyr::filter(if (isTRUE(subgroup)) {Author_pooled != "Pooled Effect" & Author_pooled != "Overall Effect" & Author_pooled != "Prediction"}
                                                  else {Author != "Pooled Effect" & Author != "Prediction"}),
          linewidth = 1, size = 1, color = color_pointinterval, limits = calc_xlim)
    }
  }
  # For common-effects: no shrinkage layer at all (tau omitted entirely)

  # Prediction layer — only when add_pred = TRUE and prediction draws exist
  if (isTRUE(add_pred)) {
    pred_data <- posterior.draws |>
      dplyr::filter(if (isTRUE(subgroup)) {
        Author_pooled == "Prediction"
      } else {
        Author == "Prediction"
      })

    if (nrow(pred_data) > 0) {
      if (pred_output == "density") {
        study.density.plot <- study.density.plot +
          ggdist::stat_slab(
            ggplot2::aes(x = x_studies, y = Author),
            data = pred_data,
            linewidth = 0.5, height = 0.9, normalize = "groups", alpha = 0.7,
            color = color_pred_outline, fill = color_pred_posterior,
            limits = calc_xlim)
      } else if (pred_output == "pointinterval") {
        study.density.plot <- study.density.plot +
          ggdist::stat_pointinterval(
            ggplot2::aes(x = x_studies, y = Author),
            data = pred_data,
            linewidth = 1, size = 1, color = color_pred_pointinterval, limits = calc_xlim)
      }
    }
  }

  return(study.density.plot)
}


#' Boxplot statistics for observed study effects (normal approximation)
#'
#' Box = interquartile range, middle = median, whiskers = 95% interval, all on
#' the model (e.g. log) scale.
#'
#' @noRd
box_stats_normal <- function(study.effects) {
  study.effects |>
    dplyr::ungroup() |>
    dplyr::mutate(se = sqrt(vi)) |>
    dplyr::transmute(
      Author,
      xmin    = yi + stats::qnorm(0.025) * se,
      xlower  = yi + stats::qnorm(0.25) * se,
      xmiddle = yi,
      xupper  = yi + stats::qnorm(0.75) * se,
      xmax    = yi + stats::qnorm(0.975) * se)
}


#' Boxplot statistics from posterior draws (model scale)
#'
#' @noRd
box_stats_draws <- function(draws, overall = FALSE) {
  draws |>
    dplyr::ungroup() |>
    dplyr::mutate(.overall = {{ overall }}) |>
    dplyr::summarise(
      xmin    = stats::quantile(theta_study, 0.025, names = FALSE),
      xlower  = stats::quantile(theta_study, 0.25, names = FALSE),
      xmiddle = stats::median(theta_study),
      xupper  = stats::quantile(theta_study, 0.75, names = FALSE),
      xmax    = stats::quantile(theta_study, 0.975, names = FALSE),
      .by = c(Author, .overall))
}


#' Move boxplot statistics to the plotting scale and squish into the x limits
#'
#' Values outside the axis limits would otherwise be dropped by the scale.
#'
#' @noRd
box_stats_scale <- function(box_df, log_scale = FALSE, xlim) {
  box_cols <- c("xmin", "xlower", "xmiddle", "xupper", "xmax")
  box_df |>
    dplyr::mutate(dplyr::across(
      dplyr::all_of(box_cols),
      \(x) pmin(pmax(if (log_scale) exp(x) else x, xlim[1]), xlim[2])))
}
