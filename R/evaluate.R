#' Evaluate fit quality
#'
#' @param signal numeric vector of calibration signal values
#' @param conc numeric vector of known concentrations
#' @param method character, fitting method name
#' @return list with method, R2, RMSE, n_points
eval_fit <- function(signal, conc, method = "linear") {
  keep <- !is.na(signal) & !is.na(conc) & is.finite(signal) & is.finite(conc)
  signal <- signal[keep]
  conc <- conc[keep]
  if (length(signal) < 2)
    stop("At least 2 valid calibration points are required for evaluation")

  fit_fn <- get_method(method)
  fit_result <- fit_fn(signal, conc)
  pred <- fit_result$predict(signal)

  ss_res <- sum((conc - pred)^2)
  ss_tot <- sum((conc - mean(conc))^2)
  r_squared <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
  rmse <- sqrt(mean((conc - pred)^2))

  list(method = tolower(trimws(method)), R2 = r_squared, RMSE = rmse, n_points = length(signal))
}

#' Compare multiple fitting methods
#'
#' @param signal numeric vector of calibration signal values
#' @param conc numeric vector of known concentrations
#' @param methods character vector of method names to compare
#' @return data.frame sorted by R2 (descending)
compare_methods <- function(signal, conc,
                            methods = c("linear", "quadratic", "log")) {
  results <- lapply(methods, function(m) {
    tryCatch(eval_fit(signal, conc, m),
             error = function(e) {
               list(method = m, R2 = NA_real_, RMSE = NA_real_,
                    n_points = length(signal), error = e$message)
             })
  })
  df <- do.call(rbind, lapply(results, function(x) {
    data.frame(method = x$method, R2 = x$R2, RMSE = x$RMSE,
               n_points = x$n_points,
               error = if (!is.null(x$error)) x$error else "",
               stringsAsFactors = FALSE)
  }))
  df[order(-df$R2), ]
}