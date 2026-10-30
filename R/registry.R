.msd_methods <- list(
  linear      = fit_linear,
  quadratic   = fit_quadratic,
  cubic       = fit_cubic,
  log         = fit_log,
  power       = fit_power,
  exponential = fit_exponential,
  "4pl"       = fit_4pl,
  "5pl"       = fit_5pl
)

#' List all available fitting methods
#' @return character vector of method names (invisibly)
list_methods <- function() {
  methods <- names(.msd_methods)
  cat("Available fitting methods:\n")
  cat("  linear      - Linear: y = a + b*x\n")
  cat("  quadratic   - Quadratic: y = a + b*x + c*x^2\n")
  cat("  cubic       - Cubic polynomial\n")
  cat("  log         - Logarithmic: y = a + b*ln(x)\n")
  cat("  power       - Power: y = a*x^b\n")
  cat("  exponential - Exponential: y = a*exp(b*x)\n")
  cat("  4pl         - Four-parameter logistic (nls-based)\n")
  cat("  5pl         - Five-parameter logistic (nls-based)\n")
  invisible(methods)
}

#' Get fitting function by name
#' @param method single fitting method name
get_method <- function(method) {
  if (!is.character(method) || length(method) != 1 || is.na(method))
    stop("method must be a single non-missing character value")
  method <- tolower(trimws(method))
  if (!method %in% names(.msd_methods)) {
    stop(sprintf("Unknown method '%s'. Available: %s",
                 method, paste(names(.msd_methods), collapse = ", ")))
  }
  .msd_methods[[method]]
}
