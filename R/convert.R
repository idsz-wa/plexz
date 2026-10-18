#' Convert signal values to concentrations with detection limits
#'
#' @param signal numeric vector of calibration signal values
#' @param conc numeric vector of known concentrations for calibration
#' @param unknown numeric vector of unknown sample signal values to convert
#' @param method character, fitting method name (see list_methods())
#' @param platform character, detection platform name (see list_platforms())
#' @return list with concentrations, status, LLOQ, ULOQ, method, platform, model_summary, fit_object
msd_to_concentration <- function(signal, conc, unknown, method = "linear",
                                 platform = "msd") {
  stopifnot(is.numeric(signal), is.numeric(conc), is.numeric(unknown))
  if (length(signal) != length(conc))
    stop("signal and conc must have the same length")

  plat <- get_platform(platform)

  keep <- !is.na(signal) & !is.na(conc) & is.finite(signal) & is.finite(conc)
  if (any(!keep)) {
    warning(sprintf("Removed %d NA or non-finite calibration point(s)", sum(!keep)))
    signal <- signal[keep]
    conc <- conc[keep]
  }
  if (length(signal) < 3)
    stop("At least 3 calibration points are required")

  if (!is.null(plat$typical_range) && length(plat$typical_range) == 2) {
    out_of_range <- sum(signal < plat$typical_range[1] | signal > plat$typical_range[2],
                        na.rm = TRUE)
    if (out_of_range > length(signal) / 2) {
      warning(sprintf(
        "Platform '%s': %d/%d calibration signals fall outside typical range [%g, %g]. Check platform assignment.",
        plat$display_name, out_of_range, length(signal),
        plat$typical_range[1], plat$typical_range[2]))
    }
  }

  fit_fn <- get_method(method)
  fit_result <- fit_fn(signal, conc)
  predicted_conc <- fit_result$predict(unknown)
  model_summary <- tryCatch(
    suppressWarnings(summary(fit_result$model)),
    error = function(e) NULL
  )

  limits <- compute_limits(signal, conc)
  status <- classify_samples(unknown, limits)

  list(
    concentrations = predicted_conc,
    status         = status,
    LLOQ           = limits$LLOQ,
    ULOQ           = limits$ULOQ,
    method         = tolower(trimws(method)),
    platform       = plat$display_name,
    model_summary  = model_summary,
    fit_object     = fit_result
  )
}

#' Platform-aware signal-to-concentration conversion with method fallback
#'
#' Convenience wrapper that auto-selects the recommended method for a platform.
#' If the data does not meet the minimum point requirement for the first
#' recommended method, it falls back through the list.
#'
#' @param signal numeric vector of calibration signal values
#' @param conc numeric vector of known concentrations for calibration
#' @param unknown numeric vector of unknown sample signal values to convert
#' @param platform character, detection platform name
#' @param method character or NULL; if NULL, auto-selects with fallback
#' @return same as msd_to_concentration()
convert_signal <- function(signal, conc, unknown, platform = "msd", method = NULL) {
  plat <- get_platform(platform)

  if (!is.null(method)) {
    return(msd_to_concentration(signal, conc, unknown, method = method, platform = platform))
  }

  n_pts <- sum(!is.na(signal) & !is.na(conc) & is.finite(signal) & is.finite(conc))

  for (m in plat$recommended_methods) {
    min_req <- if (m %in% names(.method_min_points)) .method_min_points[[m]] else 3
    if (n_pts < min_req) next
    res <- tryCatch(
      msd_to_concentration(signal, conc, unknown, method = m, platform = platform),
      error = function(e) NULL
    )
    if (!is.null(res)) {
      message(sprintf("Platform '%s': auto-selected method '%s' (%d points available)",
                      plat$display_name, m, n_pts))
      return(res)
    }
  }

  stop(sprintf(
    "Platform '%s': all recommended methods failed with %d valid points. Try specifying method manually.",
    plat$display_name, n_pts))
}

#' Multi-factor signal-to-concentration conversion with detection limits
#'
#' @param calib_data data.frame with columns: factor, signal, concentration
#' @param sample_data data.frame with columns: factor, signal
#' @param method character (applied to all factors) or named list
#'   (e.g. list("IL-6" = "4pl", "TNF-a" = "linear"))
#' @param platform character, detection platform name
#' @return data.frame with columns: factor, signal, concentration, status, LLOQ, ULOQ, method, platform
msd_multi_factor <- function(calib_data, sample_data, method = "linear",
                             platform = "msd") {
  stopifnot(is.data.frame(calib_data), is.data.frame(sample_data))

  required_calib <- c("factor", "signal", "concentration")
  required_sample <- c("factor", "signal")
  if (!all(required_calib %in% names(calib_data)))
    stop("calib_data must contain columns: ", paste(required_calib, collapse = ", "))
  if (!all(required_sample %in% names(sample_data)))
    stop("sample_data must contain columns: ", paste(required_sample, collapse = ", "))

  plat <- get_platform(platform)

  factors <- unique(sample_data$factor)
  results <- vector("list", length(factors))

  for (i in seq_along(factors)) {
    f <- factors[i]
    cal_sub <- calib_data[calib_data$factor == f, ]
    sam_sub <- sample_data[sample_data$factor == f, ]

    if (nrow(cal_sub) < 3) {
      warning(sprintf("Factor '%s' has fewer than 3 calibration points, skipping", f))
      results[[i]] <- data.frame(factor = f, signal = sam_sub$signal,
                                  concentration = NA_real_, status = "error",
                                  LLOQ = NA_real_, ULOQ = NA_real_,
                                  method = "", platform = plat$display_name,
                                  stringsAsFactors = FALSE)
      next
    }

    m <- if (is.list(method) && f %in% names(method)) method[[f]]
         else if (is.character(method)) method
         else plat$recommended_methods[1]

    res <- tryCatch(
      msd_to_concentration(cal_sub$signal, cal_sub$concentration, sam_sub$signal,
                           method = m, platform = platform),
      error = function(e) {
        warning(sprintf("Factor '%s' fit failed (%s): %s", f, m, e$message))
        NULL
      }
    )

    if (is.null(res)) {
      results[[i]] <- data.frame(factor = f, signal = sam_sub$signal,
                                  concentration = NA_real_, status = "error",
                                  LLOQ = NA_real_, ULOQ = NA_real_,
                                  method = "", platform = plat$display_name,
                                  stringsAsFactors = FALSE)
    } else {
      results[[i]] <- data.frame(factor = f, signal = sam_sub$signal,
                                  concentration = res$concentrations,
                                  status = res$status,
                                  LLOQ = res$LLOQ, ULOQ = res$ULOQ,
                                  method = res$method, platform = res$platform,
                                  stringsAsFactors = FALSE)
    }
  }

  do.call(rbind, results)
}

#' Multi-platform batch conversion
#'
#' Process multiple datasets from different platforms in one call.
#'
#' @param data_list named list of lists, each containing:
#'   calib_data, sample_data, platform, method (optional)
#' @return data.frame with all results combined, including 'batch' and 'platform' columns
multi_platform_convert <- function(data_list) {
  stopifnot(is.list(data_list))
  results <- vector("list", length(data_list))
  batch_names <- if (!is.null(names(data_list)) && all(nzchar(names(data_list))))
    names(data_list) else paste0("batch_", seq_along(data_list))

  for (i in seq_along(data_list)) {
    item <- data_list[[i]]
    if (!all(c("calib_data", "sample_data", "platform") %in% names(item))) {
      warning(sprintf("Batch '%s' missing required fields, skipping", batch_names[i]))
      results[[i]] <- data.frame(factor = NA_character_, signal = NA_real_,
                                  concentration = NA_real_, status = "error",
                                  LLOQ = NA_real_, ULOQ = NA_real_,
                                  method = "", platform = "",
                                  batch = batch_names[i], stringsAsFactors = FALSE)
      next
    }
    m <- if ("method" %in% names(item)) item$method else NULL
    p <- item$platform

    results[[i]] <- tryCatch({
      if (is.null(m)) {
        plat <- get_platform(p)
        m <- plat$recommended_methods[1]
      }
      res <- msd_multi_factor(item$calib_data, item$sample_data,
                              method = m, platform = p)
      res$batch <- batch_names[i]
      res
    }, error = function(e) {
      warning(sprintf("Batch '%s' (%s) failed: %s", batch_names[i], p, e$message))
      data.frame(factor = NA_character_, signal = NA_real_,
                  concentration = NA_real_, status = "error",
                  LLOQ = NA_real_, ULOQ = NA_real_,
                  method = "", platform = p,
                  batch = batch_names[i], stringsAsFactors = FALSE)
    })
  }

  do.call(rbind, results)
}