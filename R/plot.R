.msd_colors <- list(
  standards = "#2C73D2",
  curve     = "#D7263D",
  unknown   = "#00A878",
  limit     = "#F5A623",
  shading   = "#FBE5C8",
  grid      = "#D9D9D9"
)

#' Plot a single calibration curve with fit, unknowns, and detection limits
#'
#' @param signal numeric vector of calibration signal values
#' @param conc numeric vector of known concentrations
#' @param unknown numeric vector of unknown sample signal values (optional)
#' @param method character, fitting method name
#' @param platform character, detection platform name (affects axis label)
#' @param main plot title (defaults to method name)
#' @param log_x logical, if TRUE the x-axis is drawn on a log10 scale
#' @param show_limits logical, draw shaded out-of-range regions and LLOQ/ULOQ lines
#' @param show_grid logical, draw background grid
#' @return fit evaluation (method, R2, RMSE), invisibly
plot_calibration <- function(signal, conc, unknown = NULL, method = "linear",
                             platform = "msd", main = NULL, log_x = FALSE,
                             show_limits = TRUE, show_grid = TRUE) {
  fit_fn <- get_method(method)

  ok <- !is.na(signal) & !is.na(conc) & is.finite(signal) & is.finite(conc)
  if (any(!ok)) {
    signal <- signal[ok]
    conc <- conc[ok]
  }

  if (log_x) {
    pos <- signal > 0
    if (sum(!pos) > 0) {
      warning(sprintf("Removed %d non-positive signal point(s) for log-scale plot", sum(!pos)))
    }
    signal <- signal[pos]
    conc <- conc[pos]
    if (!is.null(unknown)) {
      unknown <- unknown[!is.na(unknown) & is.finite(unknown) & unknown > 0]
    }
  }

  if (length(signal) < 3) stop("At least 3 valid calibration points are required for plotting")

  fit_result <- fit_fn(signal, conc)
  limits <- compute_limits(signal, conc)
  pred_cal <- fit_result$predict(signal)
  ss_res <- sum((conc - pred_cal)^2)
  ss_tot <- sum((conc - mean(conc))^2)
  ev <- list(method = tolower(trimws(method)),
             R2 = if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_,
             RMSE = sqrt(mean((conc - pred_cal)^2)),
             n_points = length(signal),
             ylim = NULL, n_clipped = 0L)

  if (log_x) {
    x_cal <- log10(signal)
    y_cal <- conc
    if (!is.null(unknown) && length(unknown) > 0) {
      x_unk <- log10(unknown)
      unk_raw <- unknown
    } else {
      x_unk <- numeric(0)
      unk_raw <- numeric(0)
    }
    x_seq <- seq(min(x_cal), max(x_cal), length.out = 300)
    curve_raw <- 10^x_seq
  } else {
    x_cal <- signal
    y_cal <- conc
    if (!is.null(unknown)) {
      pos_u <- !is.na(unknown) & is.finite(unknown)
      x_unk <- unknown[pos_u]
      unk_raw <- unknown[pos_u]
    } else {
      x_unk <- numeric(0)
      unk_raw <- numeric(0)
    }
    x_seq <- seq(min(x_cal), max(x_cal), length.out = 300)
    curve_raw <- x_seq
  }

  y_seq <- fit_result$predict(curve_raw)
  y_unk <- if (length(unk_raw) > 0) fit_result$predict(unk_raw) else numeric(0)

  all_y <- c(y_cal, y_seq[is.finite(y_seq)])
  if (length(all_y) == 0) stop("No finite values available for plotting")
  ylim <- range(all_y, na.rm = TRUE)
  pad <- diff(ylim) * 0.08
  if (pad == 0) pad <- abs(ylim[1]) * 0.1 + 1
  ylim <- ylim + c(-pad, pad)

  y_unk_in <- y_unk[is.finite(y_unk) & y_unk >= ylim[1] & y_unk <= ylim[2]]
  x_unk_in <- x_unk[is.finite(y_unk) & y_unk >= ylim[1] & y_unk <= ylim[2]]
  n_clipped <- length(y_unk) - length(y_unk_in)
  ev$ylim <- ylim
  ev$n_clipped <- n_clipped

  plat <- get_platform(platform)
  sig_label <- plat$signal_label
  if (is.null(main)) main <- sprintf("%s Calibration (%s)", plat$display_name, method)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(4.5, 4.5, 3, 1.5), las = 1)

  xlab_txt <- if (log_x) sprintf("%s (%s, log scale)", sig_label, plat$display_name)
              else sprintf("%s (%s)", sig_label, plat$display_name)

  if (log_x) {
    plot(x_cal, y_cal, pch = 16, col = .msd_colors$standards, cex = 1.3,
         xlab = xlab_txt, ylab = "Concentration",
         main = main, ylim = ylim, axes = FALSE,
         xlim = range(c(x_cal, x_unk)))
    axis(2)
    ticks <- pretty(x_seq, n = 5)
    axis(1, at = ticks, labels = formatC(10^ticks, format = "g", digits = 3))
    box()
  } else {
    plot(x_cal, y_cal, pch = 16, col = .msd_colors$standards, cex = 1.3,
         xlab = xlab_txt, ylab = "Concentration",
         main = main, ylim = ylim,
         xlim = range(c(x_cal, x_unk)))
  }

  if (show_grid) grid(col = .msd_colors$grid, lty = "solid")

  if (show_limits) {
    usr <- par("usr")
    if (log_x) {
      lo_x <- if (limits$signal_lo > 0) log10(limits$signal_lo) else usr[1]
      hi_x <- if (limits$signal_hi > 0) log10(limits$signal_hi) else usr[2]
    } else {
      lo_x <- limits$signal_lo
      hi_x <- limits$signal_hi
    }
    rect(usr[1], usr[3], lo_x, usr[4],
         col = adjustcolor(.msd_colors$shading, alpha.f = 0.45), border = NA)
    rect(hi_x, usr[3], usr[2], usr[4],
         col = adjustcolor(.msd_colors$shading, alpha.f = 0.45), border = NA)
    abline(v = c(lo_x, hi_x), lty = 2, lwd = 1.5, col = .msd_colors$limit)
  }

  keep <- is.finite(y_seq)
  lines(x_seq[keep], y_seq[keep], col = .msd_colors$curve, lwd = 2)

  if (length(x_unk_in) > 0) {
    points(x_unk_in, y_unk_in, pch = 17, col = .msd_colors$unknown, cex = 1.5)
  }

  unk_label <- if (n_clipped > 0) sprintf("Unknowns (%d clipped)", n_clipped) else "Unknowns"
  legend_txt <- c("Standards",
                  sprintf("Fit (%s), R2 = %.4f", method, ev$R2),
                  unk_label, "LLOQ / ULOQ")
  legend("bottomright", legend = legend_txt,
         pch = c(16, NA, 17, NA), lty = c(NA, 1, NA, 2), lwd = c(NA, 2, NA, 1.5),
         col = c(.msd_colors$standards, .msd_colors$curve,
                 .msd_colors$unknown, .msd_colors$limit),
         bty = "n", cex = 0.85)

  invisible(ev)
}

#' Plot calibration curves for all factors as a panel grid
#'
#' @param calib_data data.frame with columns: factor, signal, concentration
#' @param sample_data data.frame with columns: factor, signal (optional)
#' @param method character or named list of methods per factor
#' @param platform character, detection platform name
#' @param ncols number of panel columns
#' @param log_x logical, apply log10 x-axis to all panels
#' @return per-factor fit evaluations, invisibly
plot_multi_factor <- function(calib_data, sample_data = NULL, method = "linear",
                              platform = "msd", ncols = 2, log_x = FALSE) {
  stopifnot(is.data.frame(calib_data))
  required <- c("factor", "signal", "concentration")
  if (!all(required %in% names(calib_data)))
    stop("calib_data must contain columns: ", paste(required, collapse = ", "))

  if (!is.null(sample_data)) {
    all_factors <- unique(c(calib_data$factor, sample_data$factor))
  } else {
    all_factors <- unique(calib_data$factor)
  }
  factors <- all_factors[all_factors %in% unique(calib_data$factor)]
  skipped <- setdiff(all_factors, factors)
  if (length(skipped) > 0) {
    warning(sprintf("Skipped %d factor(s) without calibration data: %s",
                    length(skipped), paste(skipped, collapse = ", ")))
  }
  n <- length(factors)
  if (n == 0) stop("No factors with calibration data to plot")

  ncols <- min(ncols, n)
  nrows <- ceiling(n / ncols)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(nrows, ncols), mar = c(4, 4, 2.5, 1), ask = FALSE)

  evals <- vector("list", n)
  for (i in seq_along(factors)) {
    f <- factors[i]
    cal_sub <- calib_data[calib_data$factor == f, ]
    m <- if (is.list(method) && f %in% names(method)) method[[f]]
         else if (is.character(method)) method
         else "linear"
    unk <- if (!is.null(sample_data)) sample_data$signal[sample_data$factor == f] else NULL

    evals[[i]] <- tryCatch(
      plot_calibration(cal_sub$signal, cal_sub$concentration, unknown = unk,
                       method = m, platform = platform,
                       main = sprintf("%s (%s)", f, m), log_x = log_x),
      error = function(e) {
        plot.new()
        text(0.5, 0.5, sprintf("%s: fit failed\n%s", f, e$message), cex = 0.8)
        NULL
      }
    )
  }

  invisible(evals)
}

#' Save a plot expression to an image file (PNG, PDF, JPEG, TIFF, BMP)
#'
#' @param filename output file path; format inferred from extension
#' @param expr plotting expression (evaluated inside the device)
#' @param width figure width in inches
#' @param height figure height in inches
#' @param dpi resolution in dots per inch (raster formats only)
#' @return the filename, invisibly
save_plot <- function(filename, expr, width = 8, height = 6, dpi = 300) {
  ext <- tolower(tools::file_ext(filename))
  if (ext == "") stop("filename must include an extension (e.g. .png, .pdf)")

  dev_before <- grDevices::dev.cur()

  if (ext == "pdf") {
    grDevices::pdf(filename, width = width, height = height)
  } else if (ext %in% c("png", "jpg", "jpeg", "tiff", "bmp")) {
    dev_open <- switch(ext,
      png  = function() grDevices::png(filename, width = width, height = height,
                                       units = "in", res = dpi),
      jpg  = function() grDevices::jpeg(filename, width = width, height = height,
                                        units = "in", res = dpi),
      jpeg = function() grDevices::jpeg(filename, width = width, height = height,
                                        units = "in", res = dpi),
      tiff = function() grDevices::tiff(filename, width = width, height = height,
                                        units = "in", res = dpi),
      bmp  = function() grDevices::bmp(filename, width = width, height = height,
                                       units = "in", res = dpi)
    )
    ok <- tryCatch({ dev_open(); TRUE },
                   warning = function(w) FALSE, error = function(e) FALSE)
    if (!ok) {
      while (grDevices::dev.cur() > dev_before) grDevices::dev.off()
      fallback <- sub(paste0("\\.", ext, "$"), ".pdf", filename)
      warning(sprintf("Raster device unavailable for '%s'; wrote PDF to '%s' instead",
                      filename, fallback))
      grDevices::pdf(fallback, width = width, height = height)
      filename <- fallback
    }
  } else {
    stop(sprintf("Unsupported image format '.%s'; use png, pdf, jpg, tiff, or bmp", ext))
  }

  on.exit({
    while (grDevices::dev.cur() > dev_before) grDevices::dev.off()
  }, add = TRUE)
  force(expr)
  invisible(filename)
}