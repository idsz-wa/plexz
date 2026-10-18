#' Compute detection limits (LLOQ / ULOQ) from calibration data
#'
#' @param signal numeric vector of calibration signal values
#' @param conc numeric vector of known concentrations
#' @return list with LLOQ, ULOQ (concentration values), signal_lo, signal_hi, direction
compute_limits <- function(signal, conc) {
  ord <- order(conc)
  conc_sorted <- conc[ord]
  signal_sorted <- signal[ord]
  sig_at_lloq <- signal_sorted[1]
  sig_at_uloq <- signal_sorted[length(signal_sorted)]
  direction <- if (sig_at_uloq >= sig_at_lloq) "increasing" else "decreasing"
  list(
    LLOQ      = conc_sorted[1],
    ULOQ      = conc_sorted[length(conc_sorted)],
    signal_lo = min(sig_at_lloq, sig_at_uloq),
    signal_hi = max(sig_at_lloq, sig_at_uloq),
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