#' Replace exact null-value draws with jitter so a KDE can represent them
#'
#' @description
#' Some draws carry a point mass at `null_value` alongside a continuous
#' component (e.g. RoBMA's model-averaged mu, mixing null-effect and
#' alternative-effect models). A KDE can't represent a point mass -- its
#' draws have zero variance -- so each group's exact-`null_value` draws are
#' replaced with deterministic jittered values (`stats::qnorm(stats::ppoints(.))`
#' rather than `rnorm()`, so plots stay reproducible) tightly clustered
#' around `null_value`.
#'
#' The jitter's spread is set from that *same group's own continuous draws*
#' (their range on the log scale for log-scale measures), not from the
#' plot's display range -- so the resulting spike is proportioned like a
#' feature of that group's own posterior, and groups with a wider spread
#' naturally get a wider spike. Because the null draws are only nudged, not
#' removed or resized, feeding the combined draws through one ordinary KDE
#' (`ggdist::stat_slab(normalize = "groups")`, same as every other row) makes
#' the rendered spike-to-slab proportion fall out naturally from the true
#' draw counts, with no separate height bookkeeping needed.
#'
#' @param df Data frame of draws.
#' @param group_col String naming the discrete grouping column.
#' @param x_col String naming the draws column.
#' @param null_value Numeric null value (on the same scale as `x_col`).
#' @param log_scale Logical. Whether `x_col` lives on a log-scale measure.
#' @param jitter_frac Spread of the jittered null draws, as a fraction of
#'   that group's own continuous-draws range.
#'
#' @return `df` with exact-`null_value` draws replaced by jittered values.
#'
#' @keywords internal
#' @noRd
jitter_null_draws <- function(df, group_col, x_col, null_value,
                              log_scale = FALSE, jitter_frac = 0.02) {

  df |>
    dplyr::mutate(
      "{x_col}" := {
        x <- .data[[x_col]]
        is_null <- x == null_value
        if (any(is_null) && !all(is_null)) {
          continuous <- x[!is_null]
          z <- stats::qnorm(stats::ppoints(sum(is_null)))
          if (isTRUE(log_scale)) {
            spread <- diff(range(log(continuous)))
            x[is_null] <- exp(log(null_value) + z * spread * jitter_frac)
          } else {
            spread <- diff(range(continuous))
            x[is_null] <- null_value + z * spread * jitter_frac
          }
        }
        x
      },
      .by = dplyr::all_of(group_col)
    )
}
