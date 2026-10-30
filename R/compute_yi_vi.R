#' Compute Observed Effect Sizes and Variances
#'
#' @description
#' Ensures that the data frame contains \code{yi} (observed effect size) and
#' \code{vi} (sampling variance) columns. If they already exist, the data is
#' returned unchanged. Otherwise, they are computed from canonical column names
#' based on the specified effect measure.
#'
#' This is primarily needed for one-stage models (e.g., MetaStan with raw
#' counts) where pre-calculated effect sizes are not available.
#'
#' @param data A data frame with study-level data.
#' @param measure Character string specifying the effect measure. One of
#'   "OR", "RR", "HR", "IRR", "MD", or "SMD".
#'
#' @return The input data frame with \code{yi} and \code{vi} columns added
#'   (if they were missing).
#'
#' @keywords internal
#' @noRd
compute_yi_vi <- function(data, measure = "OR") {

  # If yi and vi already exist, return as-is
  if (all(c("yi", "vi") %in% names(data))) {
    return(data)
  }

  if (measure %in% c("OR", "RR")) {
    data <- .compute_binary_yi_vi(data, measure)
  } else if (measure %in% c("MD", "SMD")) {
    data <- .compute_continuous_yi_vi(data, measure)
  } else {
    rlang::warn(
      c(
        paste0("Cannot compute yi/vi for measure '", measure, "'."),
        "i" = "Please provide yi and vi columns in your data."
      )
    )
    data$yi <- NA_real_
    data$vi <- NA_real_
  }

  data
}


#' Compute yi/vi for binary outcome measures (OR, RR)
#' @noRd
.compute_binary_yi_vi <- function(data, measure) {
  # Detect column naming convention
  if (all(c("Event_Intervention", "Event_Control",
            "N_Intervention", "N_Control") %in% names(data))) {
    r1_col <- "Event_Intervention"; r2_col <- "Event_Control"
    n1_col <- "N_Intervention";     n2_col <- "N_Control"
  } else if (all(c("r1", "r2", "n1", "n2") %in% names(data))) {
    r1_col <- "r1"; r2_col <- "r2"
    n1_col <- "n1"; n2_col <- "n2"
  } else if (all(c("rt", "rc", "nt", "nc") %in% names(data))) {
    r1_col <- "rt"; r2_col <- "rc"
    n1_col <- "nt"; n2_col <- "nc"
  } else {
    rlang::warn(
      c(
        "Cannot compute yi/vi: missing columns for binary outcome.",
        "i" = "Need event counts and sample sizes (e.g., Event_Intervention, Event_Control, N_Intervention, N_Control)."
      )
    )
    data$yi <- NA_real_
    data$vi <- NA_real_
    return(data)
  }

  # Apply 0.5 continuity correction for zero cells
  needs_correction <- data[[r1_col]] == 0 | data[[r2_col]] == 0 |
    (data[[n1_col]] - data[[r1_col]]) == 0 |
    (data[[n2_col]] - data[[r2_col]]) == 0

  r1 <- ifelse(needs_correction, data[[r1_col]] + 0.5, as.numeric(data[[r1_col]]))
  r2 <- ifelse(needs_correction, data[[r2_col]] + 0.5, as.numeric(data[[r2_col]]))
  n1 <- ifelse(needs_correction, data[[n1_col]] + 1,   as.numeric(data[[n1_col]]))
  n2 <- ifelse(needs_correction, data[[n2_col]] + 1,   as.numeric(data[[n2_col]]))

  if (measure == "OR") {
    data$yi <- log((r1 / (n1 - r1)) / (r2 / (n2 - r2)))
    data$vi <- 1/r1 + 1/(n1 - r1) + 1/r2 + 1/(n2 - r2)
  } else if (measure == "RR") {
    data$yi <- log((r1 / n1) / (r2 / n2))
    data$vi <- 1/r1 - 1/n1 + 1/r2 - 1/n2
  }

  data
}


#' Compute yi/vi for continuous outcome measures (MD, SMD)
#' @noRd
.compute_continuous_yi_vi <- function(data, measure) {
  required <- c("Mean_Intervention", "Mean_Control",
                "SD_Intervention", "SD_Control",
                "N_Intervention", "N_Control")

  if (!all(required %in% names(data))) {
    rlang::warn(
      c(
        "Cannot compute yi/vi: missing columns for continuous outcome.",
        "i" = paste("Need:", paste(required, collapse = ", "))
      )
    )
    data$yi <- NA_real_
    data$vi <- NA_real_
    return(data)
  }

  if (measure == "MD") {
    data$yi <- data$Mean_Intervention - data$Mean_Control
    data$vi <- (data$SD_Intervention^2 / data$N_Intervention) +
      (data$SD_Control^2 / data$N_Control)
  } else if (measure == "SMD") {
    # Hedges' g
    pooled_sd <- sqrt(
      ((data$N_Intervention - 1) * data$SD_Intervention^2 +
         (data$N_Control - 1) * data$SD_Control^2) /
        (data$N_Intervention + data$N_Control - 2)
    )
    j <- 1 - 3 / (4 * (data$N_Intervention + data$N_Control - 2) - 1)
    data$yi <- j * (data$Mean_Intervention - data$Mean_Control) / pooled_sd
    data$vi <- (data$N_Intervention + data$N_Control) /
      (data$N_Intervention * data$N_Control) +
      data$yi^2 / (2 * (data$N_Intervention + data$N_Control))
  }

  data
}
