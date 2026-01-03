#' Create posterior plots for Bayesian meta-analysis
#'
#' @param data The data used for the model (optional, for reference)
#' @param model A brms model object from a meta-analysis
#' @param measure Character string for effect measure ("OR", "RR", "HR", "IRR", "MD", "SMD")
#' @param study_var Character string for the study grouping variable name (default extracts from model)
#' @param incl_tau Logical, whether to include the heterogeneity (tau) plot
#' @param incl_mu_prior Logical, whether to include the prior distribution for mu
#' @param incl_tau_prior Logical, whether to include the prior distribution for tau
#' @param null_value Numeric, the null value for the effect (default from get_measure_properties)
#' @param null_range Numeric vector of length 2, range for null region (e.g., c(0.9, 1.1))
#' @param add_null_range Logical, whether to add shaded null range region
#' @param color_null_range Color for the null range shading
#' @param label_control Character string for control group label
#' @param label_intervention Character string for intervention group label
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param title_align Alignment for title ("left", "center", "right")
#' @param mu_xlim Numeric vector of length 2 for mu plot x-axis limits
#' @param tau_xlim Numeric vector of length 2 for tau plot x-axis limits
#' @param x_breaks Numeric vector for x-axis breaks on mu plot
#' @param tau_breaks Numeric vector for x-axis breaks on tau plot
#' @param color_overall_posterior Color for the overall posterior distribution
#' @param color_overall_posterior_outline Outline color for the posterior distribution
#' @param split_color_by_null Logical, whether to split posterior color by null value
#' @param color_favours_control Color for region favoring control
#' @param color_favours_intervention Color for region favoring intervention
#' @param tau_posterior_color Base color for tau posterior
#' @param tau_posterior_outline Outline color for tau posterior
#' @param tau_slab_scale Numeric, scaling factor for tau slab height relative to mu (default 0.7)
#' @param font Font family for text elements
#' @param plot_arrangement Character, "vertical" (stacked) or "horizontal" (side by side)
#'
#' @return A ggplot object (or combined plot if incl_tau = TRUE)
#' @export
#' 
#' @importFrom grDevices col2rgb rgb
#' @importFrom stats dnorm dt dcauchy
#' @importFrom ggplot2 ggplot aes theme_light element_blank element_rect element_line unit annotate geom_line geom_vline scale_fill_manual scale_x_log10 scale_x_continuous coord_cartesian labs theme margin element_text after_stat
#' @importFrom ggdist stat_slab stat_pointinterval scale_fill_ramp_discrete scale_thickness_shared
#' @importFrom patchwork plot_spacer plot_layout plot_annotation
#' @importFrom posterior as_draws_df
#' @importFrom dplyr filter group_by summarise mutate n %>%
#' @importFrom tidyr complete
#' @importFrom scales pretty_breaks
#' 
overall_plot <- function(data = NULL,
                         model,
                         measure,
                         study_var = NULL,
                         incl_tau = FALSE,
                         incl_mu_prior = FALSE,
                         incl_tau_prior = FALSE,
                         null_value = NULL,
                         null_range = NULL,
                         add_null_range = FALSE,
                         color_null_range = "#77bb41",
                         label_control = "Control",
                         label_intervention = "Intervention",
                         title = NULL,
                         subtitle = NULL,
                         title_align = "left",
                         mu_xlim = NULL,
                         tau_xlim = NULL,
                         x_breaks = NULL,
                         tau_breaks = NULL,
                         color_overall_posterior = "dodgerblue",
                         color_overall_posterior_outline = "blue",
                         split_color_by_null = FALSE,
                         color_favours_control = "firebrick",
                         color_favours_intervention = "dodgerblue",
                         tau_posterior_color = "#77bb41",
                         tau_posterior_outline = NULL,
                         tau_slab_scale = 0.7,
                         font = NULL,
                         plot_arrangement = "vertical") {
  
  # Get measure properties
  props <- get_measure_properties(measure)
  
  # Set null value from measure properties if not provided
  null_value <- null_value %||% props$null_value
  is_log_scale <- props$log_scale
  x_label <- props$x_label
  
  # Set tau x-axis label - tau represents heterogeneity (SD of random effects)
  tau_x_label <- if (measure %in% c("MD", "SMD")) {
    "Standard Deviation"
  } else {
    "Standard Deviation (Log Scale)"
  }
  
  # Set default tau outline color (darker version of base color)
  if (is.null(tau_posterior_outline)) {
    tau_rgb <- grDevices::col2rgb(tau_posterior_color)
    tau_posterior_outline <- grDevices::rgb(
      pmax(0, tau_rgb[1] * 0.5),
      pmax(0, tau_rgb[2] * 0.5),
      pmax(0, tau_rgb[3] * 0.5),
      maxColorValue = 255
    )
  }
  
  # Extract study variable name from model if not provided
  if (is.null(study_var)) {
    re_terms <- names(model$ranef)
    if (length(re_terms) > 0) {
      study_var <- re_terms[1]
    } else {
      stop("Could not automatically detect study variable.   Please specify 'study_var'.")
    }
  }
  
  # Construct the sd parameter name
  sd_param <- paste0("sd_", study_var, "__Intercept")
  
  # Extract posterior samples
  post_samples <- posterior::as_draws_df(model, c("b_Intercept", sd_param))
  
  # Extract prior information from model
  model_priors <- model$prior
  
  # --- MU PLOT ---
  
  # Transform function based on measure type
  transform_fn <- if (is_log_scale) exp else identity
  
  # Prepare mu data
  mu_samples <- post_samples$b_Intercept
  mu_transformed <- transform_fn(mu_samples)
  mu_df <- data.frame(mu = mu_transformed)
  
  # Set default x limits for mu plot with rounding
  mu_xlim <- if (!  is.null(mu_xlim)) {
    mu_xlim
  } else {
    if (is_log_scale) {
      c(0.25, 4)
    } else {
      r <- range(mu_transformed, na.rm = TRUE)
      padding <- 0.1 * diff(r)
      c(floor((r[1] - padding) * 10) / 10, ceiling((r[2] + padding) * 10) / 10)
    }
  }
  
  # Set default x breaks - include lower limit
  if (is.null(x_breaks)) {
    if (is_log_scale) {
      x_breaks <- sort(unique(c(mu_xlim[1], 0.5, 1, 2, mu_xlim[2])))
    } else {
      x_breaks <- scales::pretty_breaks(n = 5)(mu_xlim)
      if (!  mu_xlim[1] %in% x_breaks) {
        x_breaks <- sort(unique(c(mu_xlim[1], x_breaks)))
      }
    }
  } else {
    if (! mu_xlim[1] %in% x_breaks) {
      x_breaks <- sort(unique(c(mu_xlim[1], x_breaks)))
    }
  }
  
  # Determine title alignment
  title_hjust <- switch(title_align,
                        "left" = 0,
                        "center" = 0.5,
                        "right" = 1,
                        0)
  
  # Define axis line thickness
  axis_line_size <- 0.8
  
  # Base theme
  base_theme <- ggplot2::theme_light() +
    ggplot2:: theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2:: element_rect(fill = "white"),
      axis.line.x.bottom = ggplot2::element_line(color = "black", linewidth = axis_line_size),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_line(color = "black", linewidth = axis_line_size),
      axis.ticks.length = ggplot2::unit(0.15, "cm"),
      plot.title = ggplot2:: element_text(hjust = title_hjust),
      plot.subtitle = ggplot2::element_text(hjust = title_hjust),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5, "pt")
    )
  
  # Start building mu plot
  mu_plot <- ggplot2::ggplot(mu_df, ggplot2::aes(x = mu)) 
  
  # Add null range shading if requested
  if (add_null_range && !is.null(null_range)) {
    mu_plot <- mu_plot +
      ggplot2::annotate("rect", 
                        xmin = null_range[1], xmax = null_range[2],
                        ymin = -Inf, ymax = Inf,
                        fill = color_null_range, alpha = 0.2)
  }
  
  # Add prior FIRST (so it renders behind posterior) if requested
  if (incl_mu_prior) {
    # Extract intercept prior from model
    intercept_prior <- model_priors %>%
      dplyr::filter(class == "Intercept" | (class == "b" & coef == ""))
    
    if (nrow(intercept_prior) > 0) {
      prior_str <- intercept_prior$prior[1]
      mu_prior_df <- NULL
      
      # Debug: print the prior string
      message("Mu prior string:   ", prior_str)
      
      # Parse the prior string to extract distribution parameters
      if (grepl("normal", prior_str, ignore.case = TRUE) && !  grepl("student", prior_str, ignore.case = TRUE)) {
        params <- as.numeric(regmatches(prior_str, 
                                        gregexpr("-?[0-9]+\\.?[0-9]*", prior_str))[[1]])
        message("Parsed normal params: ", paste(params, collapse = ", "))
        
        if (length(params) >= 2) {
          # Create x values on the log scale for the prior
          if (is_log_scale) {
            x_seq <- seq(log(mu_xlim[1] * 0.1), log(mu_xlim[2] * 5), length.out = 500)
            prior_density <- stats::dnorm(x_seq, mean = params[1], sd = params[2])
            mu_prior_df <- data.frame(
              x = exp(x_seq),
              density = prior_density / max(prior_density) * 0.5
            )
          } else {
            x_seq <- seq(mu_xlim[1] - diff(mu_xlim) * 2, mu_xlim[2] + diff(mu_xlim) * 2, length.out = 500)
            prior_density <- stats::dnorm(x_seq, mean = params[1], sd = params[2])
            mu_prior_df <- data.frame(
              x = x_seq,
              density = prior_density / max(prior_density) * 0.5
            )
          }
        }
      } else if (grepl("student_t", prior_str, ignore.case = TRUE)) {
        params <- as.numeric(regmatches(prior_str, 
                                        gregexpr("-?[0-9]+\\.  ?[0-9]*", prior_str))[[1]])
        message("Parsed student_t params: ", paste(params, collapse = ", "))
        
        if (length(params) >= 3) {
          if (is_log_scale) {
            x_seq <- seq(log(mu_xlim[1] * 0.1), log(mu_xlim[2] * 5), length.out = 500)
            prior_density <- stats::dt((x_seq - params[2]) / params[3], df = params[1]) / params[3]
            mu_prior_df <- data.frame(
              x = exp(x_seq),
              density = prior_density / max(prior_density) * 0.5
            )
          } else {
            x_seq <- seq(mu_xlim[1] - diff(mu_xlim) * 2, mu_xlim[2] + diff(mu_xlim) * 2, length.out = 500)
            prior_density <- stats::dt((x_seq - params[2]) / params[3], df = params[1]) / params[3]
            mu_prior_df <- data.frame(
              x = x_seq,
              density = prior_density / max(prior_density) * 0.5
            )
          }
        }
      }
      
      if (!  is.null(mu_prior_df)) {
        mu_plot <- mu_plot +
          ggplot2::geom_line(data = mu_prior_df,
                             ggplot2::aes(x = x, y = density),
                             color = "grey50",
                             linewidth = 0.8,
                             linetype = "solid",
                             inherit.aes = FALSE)
      } else {
        message("mu_prior_df is NULL - prior not parsed correctly")
      }
    } else {
      message("No intercept prior found in model")
    }
  }
  
  # Add posterior distribution
  if (split_color_by_null) {
    mu_plot <- mu_plot +
      ggdist::stat_slab(ggplot2::aes(fill = ggplot2::after_stat(x < null_value)),
                        color = color_overall_posterior_outline,
                        alpha = 0.8) +
      ggplot2::scale_fill_manual(
        values = c("TRUE" = color_favours_intervention, 
                   "FALSE" = color_favours_control),
        guide = "none"
      )
  } else {
    mu_plot <- mu_plot +
      ggdist::stat_slab(fill = color_overall_posterior,
                        color = color_overall_posterior_outline,
                        alpha = 0.8)
  }
  
  # Add point interval
  mu_plot <- mu_plot +
    ggdist:: stat_pointinterval(.width = c(.66, .80, .95), 
                                color = "black",
                                point_size = 2)
  
  # Add null reference line (black, solid)
  mu_plot <- mu_plot +
    ggplot2::geom_vline(xintercept = null_value, linetype = "solid", color = "black", linewidth = 0.8)
  
  # Add direction labels (same y-level, bold) - at top
  label_y <- 0.95
  if (is_log_scale) {
    left_x <- mu_xlim[1] * 1.4
    right_x <- mu_xlim[2] * 0.65
  } else {
    range_width <- diff(mu_xlim)
    left_x <- mu_xlim[1] + range_width * 0.12
    right_x <- mu_xlim[2] - range_width * 0.12
  }
  
  mu_plot <- mu_plot +
    ggplot2::annotate("text", x = left_x, y = label_y, 
                      label = paste0("Favours\n", label_intervention),
                      fontface = "bold.italic", size = 3, color = "grey30", vjust = 1) +
    ggplot2:: annotate("text", x = right_x, y = label_y, 
                       label = paste0("Favours\n", label_control),
                       fontface = "bold.italic", size = 3, color = "grey30", vjust = 1) +
    ggplot2:: annotate("segment", x = mu_xlim[1], xend = mu_xlim[2], y = Inf, yend = Inf,
                       linewidth = axis_line_size, color = "grey60"
    )
  
  # Scale and coordinate system - minimal gap at bottom (0.03)
  if (is_log_scale) {
    mu_plot <- mu_plot +
      ggplot2::scale_x_log10(breaks = x_breaks, 
                             limits = mu_xlim,
                             labels = function(x) sprintf("%.2g", x),
                             expand = c(0, 0)) +
      ggplot2:: coord_cartesian(xlim = mu_xlim, ylim = c(0.03, 1), clip = "off")
  } else {
    mu_plot <- mu_plot +
      ggplot2::scale_x_continuous(breaks = x_breaks, 
                                  limits = mu_xlim,
                                  labels = function(x) sprintf("%.2g", x),
                                  expand = c(0, 0)) +
      ggplot2::coord_cartesian(xlim = mu_xlim, ylim = c(0.03, 1), clip = "off")
  }
  
  mu_plot <- mu_plot +
    ggdist::scale_thickness_shared() +
    base_theme +
    ggplot2::labs(x = x_label, 
                  y = NULL,
                  title = expression(paste("Overall Effect (", mu, ")")))
  
  # Apply custom font if specified
  if (!is.null(font)) {
    mu_plot <- mu_plot +
      ggplot2::theme(text = ggplot2::element_text(family = font))
  }
  
  # --- TAU PLOT (if requested) ---
  
  if (incl_tau) {
    # Extract tau samples
    tau_samples <- post_samples[[sd_param]]
    tau_df <- data.frame(tau = tau_samples)
    
    # Define heterogeneity cutpoints
    tau_cuts <- c(0, 0.1, 0.5, 1, Inf)
    tau_labels <- c("Low", "Reasonable", "Fairly high", "Fairly extreme")
    
    # Categorize tau samples and calculate percentages
    tau_df$category <- cut(tau_df$tau, 
                           breaks = tau_cuts, 
                           labels = tau_labels,
                           include.lowest = TRUE,
                           right = FALSE)
    
    # Calculate percentage in each category
    tau_percentages <- tau_df %>%
      dplyr::group_by(category) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(pct = round(n / sum(n) * 100, 1)) %>%
      tidyr::complete(category = factor(tau_labels, levels = tau_labels), 
                      fill = list(n = 0, pct = 0))
    
    # Create labels with percentages
    pct_low <- tau_percentages$pct[tau_percentages$category == "Low"]
    pct_reasonable <- tau_percentages$pct[tau_percentages$category == "Reasonable"]
    pct_fairly_high <- tau_percentages$pct[tau_percentages$category == "Fairly high"]
    pct_extreme <- tau_percentages$pct[tau_percentages$category == "Fairly extreme"]
    
    # Set default tau x limits with rounding
    tau_xlim <- if (!  is.null(tau_xlim)) {
      tau_xlim
    } else {
      r <- range(tau_df$tau, na.rm = TRUE)
      padding <- 0.1 * diff(r)
      c(max(0, floor((r[1] - padding) * 10) / 10), ceiling((r[2] + padding) * 10) / 10)
    }
    
    # Set default tau breaks - include lower limit and heterogeneity cutpoints
    if (is.null(tau_breaks)) {
      # Base heterogeneity reference points
      reference_points <- c(0, 0.10, 0.25, 0.50, 0.75, 1.00)
      
      # Add intermediate breaks between 1 and tau_xlim[2] if needed
      if (tau_xlim[2] > 1) {
        # Add breaks at 0.5 intervals above 1 (1.5, 2.0, 2.5, etc.)
        extra_breaks <- seq(1.5, tau_xlim[2], by = 0.5)
        extra_breaks <- extra_breaks[extra_breaks < tau_xlim[2]]  # Don't duplicate upper limit
        reference_points <- c(reference_points, extra_breaks)
      }
      
      tau_breaks <- reference_points[reference_points >= tau_xlim[1] & reference_points <= tau_xlim[2]]
      
      # Always include the limits
      tau_breaks <- sort(unique(c(tau_xlim[1], tau_breaks, tau_xlim[2])))
    } else {
      # User provided breaks - ensure limits are included
      if (! tau_xlim[1] %in% tau_breaks) {
        tau_breaks <- sort(unique(c(tau_xlim[1], tau_breaks)))
      }
      if (! tau_xlim[2] %in% tau_breaks) {
        tau_breaks <- sort(unique(c(tau_breaks, tau_xlim[2])))
      }
    }
    
    # Build tau plot with scaled slab height
    tau_plot <- ggplot2::ggplot(tau_df, ggplot2:: aes(x = tau)) +
      ggdist::stat_slab(
        ggplot2::aes(fill_ramp = ggplot2::after_stat(cut(x, 
                                                         breaks = c(-Inf, 0.1, 0.5, 1, Inf),
                                                         labels = c("Low", "Reasonable", "Fairly high", "Fairly extreme")))),
        fill = tau_posterior_color,
        color = tau_posterior_outline,
        linewidth = 0.8,
        scale = tau_slab_scale  # Scale the slab height
      ) +
      ggdist::scale_fill_ramp_discrete(
        range = c(0.8, 0.2),
        guide = "none"
      ) +
      ggdist::stat_pointinterval(.width = c(.66, .80, .95), 
                                 color = "black",
                                 point_size = 2)
    
    # Add tau prior if requested
    if (incl_tau_prior) {
      sd_prior <- model_priors %>%
        dplyr::filter(class == "sd")
      
      if (nrow(sd_prior) > 0) {
        prior_str <- sd_prior$prior[1]
        tau_prior_df <- NULL
        
        message("Tau prior string:  ", prior_str)
        
        if (grepl("cauchy|half_cauchy", prior_str, ignore.case = TRUE)) {
          params <- as.numeric(regmatches(prior_str, 
                                          gregexpr("-? [0-9]+\\. ?[0-9]*", prior_str))[[1]])
          if (length(params) >= 2) {
            x_seq <- seq(0.001, tau_xlim[2] * 1.5, length.out = 500)
            prior_density <- stats::dcauchy(x_seq, location = params[1], scale = params[2])
            prior_density <- prior_density * 2
            tau_prior_df <- data.frame(
              x = x_seq,
              density = prior_density / max(prior_density) * tau_slab_scale * 0.5
            )
          } else if (length(params) >= 1) {
            x_seq <- seq(0.001, tau_xlim[2] * 1.5, length.out = 500)
            prior_density <- stats::dcauchy(x_seq, location = 0, scale = params[1])
            prior_density <- prior_density * 2
            tau_prior_df <- data.frame(
              x = x_seq,
              density = prior_density / max(prior_density) * tau_slab_scale * 0.5
            )
          }
        } else if (grepl("normal|half_normal", prior_str, ignore.case = TRUE)) {
          params <- as.numeric(regmatches(prior_str, 
                                          gregexpr("-?  [0-9]+\\.?  [0-9]*", prior_str))[[1]])
          if (length(params) >= 2) {
            x_seq <- seq(0.001, tau_xlim[2] * 1.5, length.out = 500)
            prior_density <- stats::dnorm(x_seq, mean = params[1], sd = params[2])
            prior_density <- prior_density * 2
            tau_prior_df <- data.frame(
              x = x_seq,
              density = prior_density / max(prior_density) * tau_slab_scale * 0.5
            )
          }
        } else if (grepl("student_t", prior_str, ignore.case = TRUE)) {
          params <- as.numeric(regmatches(prior_str, 
                                          gregexpr("-?[0-9]+\\.?[0-9]*", prior_str))[[1]])
          if (length(params) >= 3) {
            x_seq <- seq(0.001, tau_xlim[2] * 1.5, length.out = 500)
            prior_density <- stats::dt((x_seq - params[2]) / params[3], df = params[1]) / params[3]
            prior_density <- prior_density * 2
            tau_prior_df <- data.frame(
              x = x_seq,
              density = prior_density / max(prior_density) * tau_slab_scale * 0.5
            )
          }
        }
        
        if (!  is.null(tau_prior_df)) {
          tau_plot <- tau_plot +
            ggplot2::geom_line(data = tau_prior_df,
                               ggplot2::aes(x = x, y = density),
                               color = "grey50",
                               linewidth = 0.8,
                               linetype = "solid",
                               inherit.aes = FALSE)
        }
      }
    }
    
    # Add heterogeneity category reference lines (dashed, grey)
    tau_plot <- tau_plot +
      ggplot2::geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
      ggplot2:: geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40", linewidth = 0.5) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.5)
    
    # Calculate center position for "Very High" label (between 1 and tau_xlim[2])
    very_high_x <- (1 + tau_xlim[2]) / 2
    
    # Add labels with percentages (same y-level, bold) - at top
    # Only include labels if their region falls within tau_xlim
    
    # Low:   region is 0 to 0.1, label at 0.05
    # Show if tau_xlim[1] < 0.1 (left boundary doesn't cut off entire region) AND tau_xlim[2] > 0 (right boundary includes start)
    if (tau_xlim[1] < 0.1 && tau_xlim[2] > 0) {
      tau_plot <- tau_plot +
        ggplot2::annotate("text", x = 0.05, y = 0.95, 
                          label = paste0("Low\n(", pct_low, "%)"), 
                          fontface = "bold.italic", size = 2.6, color = "grey30", vjust = 1)
    }
    
    # Moderate: region is 0.1 to 0.5, label at 0.3
    # Show if tau_xlim[1] < 0.5 (left boundary doesn't cut off entire region) AND tau_xlim[2] > 0.1 (right boundary includes start)
    if (tau_xlim[1] < 0.5 && tau_xlim[2] > 0.1) {
      tau_plot <- tau_plot +
        ggplot2:: annotate("text", x = 0.3, y = 0.95, 
                           label = paste0("Moderate\n(", pct_reasonable, "%)"), 
                           fontface = "bold.italic", size = 2.6, vjust = 1)
    }
    
    # High: region is 0.5 to 1.0, label at 0.75
    # Show if tau_xlim[1] < 1.0 (left boundary doesn't cut off entire region) AND tau_xlim[2] > 0.5 (right boundary includes start)
    if (tau_xlim[1] < 1.0 && tau_xlim[2] > 0.5) {
      tau_plot <- tau_plot +
        ggplot2::annotate("text", x = 0.75, y = 0.95, 
                          label = paste0("High\n(", pct_fairly_high, "%)"), 
                          fontface = "bold.italic", size = 2.6, vjust = 1)
    }
    
    # Very High: region is 1.0+, label at very_high_x
    # Show if tau_xlim[2] > 1.0 (right boundary extends past 1.0)
    if (tau_xlim[2] > 1.0) {
      tau_plot <- tau_plot +
        ggplot2::annotate("text", x = very_high_x, y = 0.95, 
                          label = paste0("Very High\n(", pct_extreme, "%)"), 
                          fontface = "bold.italic", size = 2.6, vjust = 1)
    }
    
    tau_plot <- tau_plot +
      ggplot2::annotate("segment", x = tau_xlim[1], xend = tau_xlim[2], y = Inf, yend = Inf,
                        linewidth = axis_line_size, color = "grey60"
      )
    
    # Scale and theme for tau plot - minimal gap at bottom
    tau_plot <- tau_plot +
      ggplot2:: scale_x_continuous(breaks = tau_breaks, 
                                   limits = tau_xlim,
                                   labels = function(x) sprintf("%.2g", x),
                                   expand = c(0, 0)) +
      ggplot2::coord_cartesian(xlim = tau_xlim, ylim = c(0.03, 1), clip = "off") +
      ggdist::scale_thickness_shared() +
      base_theme +
      ggplot2::labs(x = tau_x_label, 
                    y = NULL,
                    title = expression(paste("Heterogeneity (", tau, ")")))
    
    # Apply custom font if specified
    if (!is.null(font)) {
      tau_plot <- tau_plot +
        ggplot2::theme(text = ggplot2::element_text(family = font))
    }
    
    # Combine plots using patchwork
    if (plot_arrangement == "horizontal") {
      combined_plot <- mu_plot + patchwork::plot_spacer() + tau_plot +
        patchwork::plot_layout(ncol = 3, widths = c(1, 0.05, 1))
    } else {
      combined_plot <- mu_plot / tau_plot +
        patchwork::plot_layout(ncol = 1)
    }
    
    # Add title and subtitle using patchwork
    if (!  is.null(title) || ! is.null(subtitle)) {
      annotation_theme <- ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = title_hjust, face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(hjust = title_hjust, size = 11)
      )
      
      if (! is.null(font)) {
        annotation_theme <- annotation_theme +
          ggplot2::theme(
            text = ggplot2::element_text(family = font),
            plot.title = ggplot2::element_text(hjust = title_hjust, face = "bold", size = 14, family = font),
            plot.subtitle = ggplot2::element_text(hjust = title_hjust, size = 11, family = font)
          )
      }
      
      combined_plot <- combined_plot +
        patchwork::plot_annotation(
          title = title,
          subtitle = subtitle,
          theme = annotation_theme
        )
    }
    
    return(combined_plot)
    
  } else {
    # Return just the mu plot with title/subtitle if provided
    if (!is.null(title) || !is.null(subtitle)) {
      mu_plot <- mu_plot +
        ggplot2::labs(title = if (!is.null(title)) title else expression(paste("Overall Effect (", mu, ")")),
                      subtitle = subtitle)
    }
    return(mu_plot)
  }
}