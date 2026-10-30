#' Compute detection limits (LLOQ / ULOQ) from calibration data
#'
#' @param signal numeric vector of calibration signal values
#' @param conc numeric vector of known concentrations
#' @return list with LLOQ, ULOQ (concentration values), signal_lo, signal_hi, direction
compute_limits <- function(signal, conc) {
  if (!is.numeric(signal) || !is.numeric(conc) || length(signal) != length(conc))
    stop("signal and conc must be numeric vectors of the same length")
  keep <- is.finite(signal) & is.finite(conc)
  if (sum(keep) < 2)
    stop("At least 2 finite calibration points are required")
  signal <- signal[keep]
  conc <- conc[keep]
  if (length(unique(conc)) < 2)
    stop("At least 2 distinct concentration values are required")

  ord <- order(conc)
  conc_sorted <- conc[ord]
  signal_sorted <- signal[ord]
  group_signal <- vapply(split(signal_sorted, conc_sorted), median, numeric(1))
  sig_at_lloq <- group_signal[1]
  sig_at_uloq <- group_signal[length(group_signal)]
  direction <- if (sig_at_uloq >= sig_at_lloq) "increasing" else "decreasing"

  diffs <- diff(group_signal)
  if (any(diffs != 0) && any(sign(diffs[diffs != 0]) !=
                             if (direction == "increasing") 1 else -1)) {
    warning("Calibration signals are not monotonic with concentration; using the observed signal range")
  }

  list(
    LLOQ      = conc_sorted[1],
    ULOQ      = conc_sorted[length(conc_sorted)],
    signal_lo = min(signal_sorted),
    signal_hi = max(signal_sorted),
    direction = direction
  )
}

#' Classify samples relative to detection limits
#'
#' @param unknown_signal numeric vector of sample signal values
#' @param limits list from compute_limits()
#' @return character vector: "in_range", "below_LLOQ", "above_ULOQ", or "unknown" for non-finite signals
classify_samples <- function(unknown_signal, limits) {
  valid <- !is.na(unknown_signal) & is.finite(unknown_signal)
  status <- rep("unknown", length(unknown_signal))

  if (limits$direction == "increasing") {
    status[valid & unknown_signal < limits$signal_lo] <- "below_LLOQ"
    status[valid & unknown_signal > limits$signal_hi] <- "above_ULOQ"
    status[valid & unknown_signal >= limits$signal_lo & unknown_signal <= limits$signal_hi] <- "in_range"
  } else {
    status[valid & unknown_signal > limits$signal_hi] <- "below_LLOQ"
    status[valid & unknown_signal < limits$signal_lo] <- "above_ULOQ"
    status[valid & unknown_signal >= limits$signal_lo & unknown_signal <= limits$signal_hi] <- "in_range"
  }
  status
}

#' Summarize detection limits and sample status per factor
#'
#' @param result data.frame returned by msd_multi_factor()
#' @return data.frame with one row per factor: LLOQ, ULOQ, n_in_range, n_below, n_above, n_error
summarize_limits <- function(result) {
  stopifnot(is.data.frame(result))
  required <- c("factor", "status", "LLOQ", "ULOQ")
  if (!all(required %in% names(result)))
    stop("result must contain columns: ", paste(required, collapse = ", "))

  factors <- unique(result$factor)
  rows <- lapply(factors, function(f) {
    sub <- result[result$factor == f, ]
    data.frame(
      factor     = f,
      LLOQ       = sub$LLOQ[1],
      ULOQ       = sub$ULOQ[1],
      n_total    = nrow(sub),
      n_in_range = sum(sub$status == "in_range", na.rm = TRUE),
      n_below    = sum(sub$status == "below_LLOQ", na.rm = TRUE),
      n_above    = sum(sub$status == "above_ULOQ", na.rm = TRUE),
      n_unknown  = sum(sub$status == "unknown" | is.na(sub$status)),
      n_error    = sum(sub$status == "error", na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
