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
#' @param study_output "density", "boxplot", or "sq_whisker" for the
#'   per-study rows. Boxplots and square-whiskers are horizontal: the
#'   interval spans the 95% interval; boxplots also show the interquartile
#'   range as a box.
#' @param pooled_output "density", "boxplot", "sq_whisker", or "diamond" for
#'   the Pooled/Overall Effect row(s) -- independent of `study_output`.
#' @param shrinkage_output "density" or "pointinterval".
#' @param pred_output "density", "pointinterval", "boxplot", "sq_whisker", or
#'   "diamond".
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
                                   color_shrinkage_pointinterval = "purple",
                                   color_shrinkage_fill = NA,
                                   split_color_by_null = FALSE,
                                   color_favours_control = "firebrick",
                                   color_favours_intervention = "dodgerblue",
                                   label_control = "Control",
                                   label_intervention = "Intervention",
                                   add_arm_labels = TRUE,
                                   reverse_arm_labels = FALSE,
                                   study_output = "density",
                                   pooled_output = "density",
                                   pooled_labels = "Pooled Effect",
                                   shrinkage_output = "density",
                                   add_pred = FALSE,
                                   add_pred_subgroup = FALSE,
                                   pred_output = "density",
                                   pred_labels = "Prediction",
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
                                   section_break_labels = character(0),
                                   font = NULL){

  box_styles <- c("boxplot", "sq_whisker", "diamond")

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
    color_shrinkage_pointinterval <- palette_colors[6]
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

  # Box-style glyph layer builder (boxplot / sq_whisker / diamond)
  box_glyph_layer <- function(box_df, style, outline_color, width = 0.35) {
    mapping <- if (isTRUE(split_color_by_null)) {
      ggplot2::aes(y = Author, x = xmiddle, xmin = xmin, xlower = xlower,
                   xmiddle = xmiddle, xupper = xupper, xmax = xmax, fill = .above_null)
    } else {
      ggplot2::aes(y = Author, x = xmiddle, xmin = xmin, xlower = xlower,
                   xmiddle = xmiddle, xupper = xupper, xmax = xmax, fill = .fill)
    }
    position <- ggplot2::position_nudge(y = 0.5)
    switch(style,
      boxplot = ggplot2::geom_boxplot(
        mapping, data = box_df, stat = "identity", orientation = "y", width = width,
        position = position, alpha = 0.7,
        colour = outline_color, linewidth = 0.5, inherit.aes = FALSE),
      sq_whisker = ggplot2::geom_pointrange(
        utils::modifyList(mapping, ggplot2::aes(colour = ggplot2::after_scale(fill))),
        data = box_df, orientation = "y",
        position = position, shape = 22, size = 0.7,
        linewidth = 0.8, inherit.aes = FALSE),
      diamond = geom_forest_diamond(
        mapping, data = box_df, height = width * 0.7,
        position = position,
        colour = outline_color, linewidth = 0.5, inherit.aes = FALSE))
  }

  # Study-row box_df (boxplot / sq_whisker)
  study_box_layer <- if (study_output %in% box_styles) {
    study_box_df <- study.effects |>
      dplyr::filter(!is.na(yi), !is.na(vi)) |>
      box_stats_normal() |>
      box_stats_scale(log_scale = isTRUE(props$log_scale), xlim = calc_xlim) |>
      dplyr::mutate(
        .fill = color_study_posterior,
        .above_null = xmiddle > null_value)
    box_glyph_layer(study_box_df, study_output, color_study_posterior_outline)
  }

  # Pooled-row box_df (boxplot / sq_whisker / diamond)
  pooled_box_layer <- if (pooled_output %in% box_styles) {
    pooled_rows <- if (isTRUE(subgroup)) {
      posterior.draws |> dplyr::filter(Author_pooled %in% c("Pooled Effect", "Overall Effect"))
    } else {
      posterior.draws |> dplyr::filter(Author %in% pooled_labels)
    }
    pooled_box_df <- pooled_rows |>
      box_stats_draws(overall = if (isTRUE(subgroup)) Author_pooled == "Overall Effect" else FALSE) |>
      box_stats_scale(log_scale = isTRUE(props$log_scale), xlim = calc_xlim) |>
      dplyr::mutate(
        .fill = dplyr::if_else(.overall, color_overall_posterior, color_pooled_posterior),
        .above_null = xmiddle > null_value) |>
      dplyr::select(-.overall)
    box_glyph_layer(pooled_box_df, pooled_output, color_study_posterior_outline)
  }

  box_ylim <- if (study_output != "density" && pooled_output != "density") {
    c(1, dplyr::n_distinct(posterior.draws$Author) + 0.9)
  }

  pooled_ref_lines <- if (isTRUE(props$log_scale)) {
    list(
      ggplot2::geom_vline(xintercept = exp(pooled_effect["estimate"]), color = "grey60", linewidth = 1),
      ggplot2::geom_vline(xintercept = exp(pooled_effect[c("lower", "upper")]), color = "grey60", linetype = 2))
  } else {
    list(
      ggplot2::geom_vline(xintercept = pooled_effect["estimate"], color = "grey60", linewidth = 1),
      ggplot2::geom_vline(xintercept = pooled_effect[c("lower", "upper")], color = "grey60", linetype = 2))
  }

  study.density.plot <-
    ggplot2::ggplot(ggplot2::aes(y = Author), data = posterior.draws) +
    # Add ROPE shading first
    {if (isTRUE(add_rope)) {
      ggplot2::annotate("rect",
                        xmin = rope_range[1], xmax = rope_range[2],
                        ymin = -Inf, ymax = Inf,
                        fill = scales::alpha(rope_color, 0.3),
                        color = NA)
    }} +
    # Median/CrI reference lines next
    pooled_ref_lines +
    {if (study_output != "density") NULL else if (isTRUE(split_color_by_null)) {
      list(
        ggdist::stat_slab(
          ggplot2::aes(xdist = xdist, fill = ggplot2::after_stat(x > null_value)),
          slab_linewidth = 0.5, alpha = 0.7, limits = calc_xlim, height = 0.9,
          normalize = "groups",
          data = study.effects, colour = color_study_posterior_outline))
    } else {
      ggdist::stat_slab(
        ggplot2::aes(xdist = xdist),
        slab_linewidth = 0.5, alpha = 0.7, limits = calc_xlim, height = 0.9,
        normalize = "groups",
        data = study.effects, colour = color_study_posterior_outline,
        fill = color_study_posterior)
    }} +
    {if (pooled_output != "density") NULL else if (isTRUE(split_color_by_null)) {
      list(
        ggdist::stat_slab(
          ggplot2::aes(
            x = x_studies, y = Author, fill = ggplot2::after_stat(x > null_value)),
          data = posterior.draws |>
            dplyr::filter(if (isTRUE(subgroup)) {
              Author_pooled == "Pooled Effect"
            } else {
              Author %in% pooled_labels
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
            Author %in% pooled_labels
          }),
        fill = color_pooled_posterior, height = 0.9, normalize = "groups",
        limits = calc_xlim)
    }} +
    # Add overall effect slab when subgroup = TRUE
    {if (isTRUE(subgroup) && pooled_output == "density") {
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
    study_box_layer +
    pooled_box_layer +
    {if (isTRUE(split_color_by_null)) {
      ggplot2::scale_fill_manual(
        values = c(
          "FALSE" = color_favours_intervention,
          "TRUE"  = color_favours_control),
        guide = "none")
    } else if (study_output %in% box_styles || pooled_output %in% box_styles) {
      ggplot2::scale_fill_identity()
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
          xmin = calc_xlim[1] - 0.01, xmax = calc_xlim[2], ymin = -Inf, ymax = Inf))
    }} +
    ggplot2::labs(x = x_label %||% props$x_label) +
    ggplot2::ylab(NULL) +
    ggplot2::geom_hline(yintercept = 1, color = "black", linewidth = 0.75) +
    {if (length(section_break_labels) > 0) {
      author_levels_all <- levels(posterior.draws$Author)
      level_index <- match(section_break_labels, author_levels_all)
      break_positions <- length(author_levels_all) - level_index + 2
      purrr::map(break_positions[!is.na(break_positions)], ~ {
        ggplot2::geom_hline(yintercept = .x, color = "black", linewidth = 0.4)
      })
    }}

  if (isTRUE(props$log_scale)){
    study.density.plot <- study.density.plot +
      ggplot2::scale_x_log10(breaks = breaks, expand = c(0, 0), limits = calc_xlim)
  } else {
    study.density.plot <- study.density.plot +
      ggplot2::scale_x_continuous(breaks = breaks, expand = c(0, 0), limits = calc_xlim)
  }

  if(isFALSE(subgroup)){
    study.density.plot <- study.density.plot +
      ggplot2::scale_y_discrete(expand = c(0, 0), limits = rev(levels(posterior.draws$Author)))
  } else if (isTRUE(subgroup)){
    # Fix the rows explicitly so empty subgroup header rows are kept
    study.density.plot <- study.density.plot +
      ggplot2::scale_y_discrete(expand = c(0, 0), limits = levels(posterior.draws$Author))
  }

  if(is.null(shrinkage_output))"density" else shrinkage_output

  # Shrinkage layer
  if (effect_type == "random_effect") {
    if (shrinkage_output == "density") {
      study.density.plot <- study.density.plot +
        ggdist::stat_slab(
          ggplot2::aes(x = x_studies, y = Author),
          data = posterior.draws |> dplyr::filter(if (isTRUE(subgroup)) {Author_pooled != "Pooled Effect" & Author_pooled != "Overall Effect" & Author_pooled != "Prediction"}
                                                  else {!Author %in% c(pooled_labels, pred_labels)}),
          linewidth = 0.5, scale = 0.6, normalize = "groups",
          color = color_shrinkage_outline, fill = color_shrinkage_fill, limits = calc_xlim)
    } else if (shrinkage_output == "pointinterval") {
      study.density.plot <- study.density.plot +
        ggdist::stat_pointinterval(
          ggplot2::aes(x = x_studies, y = Author),
          data = posterior.draws |> dplyr::filter(if (isTRUE(subgroup)) {Author_pooled != "Pooled Effect" & Author_pooled != "Overall Effect" & Author_pooled != "Prediction"}
                                                  else {!Author %in% c(pooled_labels, pred_labels)}),
          linewidth = 1, size = 1, color = color_shrinkage_pointinterval, limits = calc_xlim)
    }
  }

  # Prediction layer
  if (isTRUE(add_pred) || !identical(pred_labels, "Prediction")) {
    pred_data <- posterior.draws |>
      dplyr::filter(if (isTRUE(subgroup)) {
        Author_pooled == "Prediction"
      } else {
        Author %in% pred_labels
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
      } else if (pred_output %in% box_styles) {
        pred_box_df <- pred_data |>
          box_stats_draws() |>
          box_stats_scale(log_scale = isTRUE(props$log_scale), xlim = calc_xlim)

        mapping <- ggplot2::aes(y = Author, x = xmiddle, xmin = xmin, xlower = xlower,
                                xmiddle = xmiddle, xupper = xupper, xmax = xmax)
        pred_position <- ggplot2::position_nudge(y = 0.5)
        study.density.plot <- study.density.plot +
          switch(pred_output,
            boxplot = ggplot2::geom_boxplot(
              mapping, data = pred_box_df, stat = "identity", orientation = "y", width = 0.35,
              position = pred_position, alpha = 0.7,
              fill = color_pred_posterior, colour = color_pred_outline, linewidth = 0.5,
              inherit.aes = FALSE),
            sq_whisker = ggplot2::geom_pointrange(
              mapping, data = pred_box_df, orientation = "y",
              position = pred_position, shape = 22, size = 0.7,
              fill = color_pred_posterior, colour = color_pred_posterior, linewidth = 0.8,
              inherit.aes = FALSE),
            diamond = geom_forest_diamond(
              mapping, data = pred_box_df, height = 0.35 * 0.7,
              position = pred_position,
              fill = color_pred_posterior, colour = color_pred_outline, linewidth = 0.5,
              inherit.aes = FALSE))
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
      xmin = yi + stats::qnorm(0.025) * se,
      xlower = yi + stats::qnorm(0.25) * se,
      xmiddle = yi,
      xupper = yi + stats::qnorm(0.75) * se,
      xmax = yi + stats::qnorm(0.975) * se)
}


#' Boxplot statistics from posterior draws
#'
#' @noRd
box_stats_draws <- function(draws, overall = FALSE) {
  draws |>
    dplyr::ungroup() |>
    dplyr::mutate(.overall = {{ overall }}) |>
    dplyr::summarise(
      xmin = stats::quantile(theta_study, 0.025, names = FALSE),
      xlower = stats::quantile(theta_study, 0.25, names = FALSE),
      xmiddle = stats::median(theta_study),
      xupper= stats::quantile(theta_study, 0.75, names = FALSE),
      xmax = stats::quantile(theta_study, 0.975, names = FALSE),
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


# Diamond geom

#' @noRd
GeomForestDiamond <- ggplot2::ggproto("GeomForestDiamond", ggplot2::Geom,
  required_aes = c("x", "xmin", "xmax", "y"),
  default_aes = ggplot2::aes(
    colour = "black", fill = "grey20", linewidth = 0.5, linetype = 1,
    height = 0.35, alpha = NA),
  draw_key = ggplot2::draw_key_polygon,
  setup_data = function(data, params) {
    data$height <- data$height %||% params$height %||% 0.35
    data
  },
  draw_panel = function(data, panel_params, coord, lineend = "butt") {
    vertices <- purrr::map_dfr(seq_len(nrow(data)), function(i) {
      row <- data[i, ]
      tibble::tibble(
        x = c(row$xmin, row$x, row$xmax, row$x),
        y = c(row$y, row$y + row$height, row$y, row$y - row$height),
        group = i,
        colour = row$colour, fill = row$fill, linewidth = row$linewidth,
        linetype = row$linetype, alpha = row$alpha)
    })
    ggplot2::GeomPolygon$draw_panel(vertices, panel_params, coord, lineend = lineend)
  }
)

#' Meta-analysis diamond glyph
#'
#' A horizontal diamond spanning `xmin` to `xmax`, pointed at `x`, placed at
#' `y`. Used for pooled-effect and prediction rows in [bayes_forest()].
#'
#' @noRd
geom_forest_diamond <- function(mapping = NULL, data = NULL, stat = "identity",
                                position = "identity", ..., na.rm = FALSE,
                                show.legend = NA, inherit.aes = TRUE) {
  ggplot2::layer(
    data = data, mapping = mapping, stat = stat, geom = GeomForestDiamond,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = rlang::list2(na.rm = na.rm, ...))
}
