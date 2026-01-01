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

#' Extract Intercept Draws from Model
#'
#' @description
#' Extracts posterior draws for the intercept parameter and transforms them
#' according to the effect measure.
#'
#' @param model A fitted brmsfit object
#' @param measure Character string. Effect measure ("OR", "RR", "HR", "IRR", "MD", "SMD")
#'
#' @return Numeric vector of transformed posterior draws
#'
#' @keywords internal
#' @noRd
extract_intercept_draws <- function(model, measure) {
  draws <- tidybayes::spread_draws(model, b_Intercept)
  
  if (measure %in% c("OR", "RR", "HR", "IRR")) {
    draws$x <- exp(draws$b_Intercept)
  } else {
    draws$x <- draws$b_Intercept
  }
  
  draws$x
}