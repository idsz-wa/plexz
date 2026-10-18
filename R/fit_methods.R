#' Linear fit: y = a + b*x
fit_linear <- function(signal, conc) {
  stopifnot(is.numeric(signal), is.numeric(conc), length(signal) == length(conc))
  model <- lm(conc ~ signal)
  list(model = model, name = "linear",
       predict = function(s) as.numeric(predict(model, newdata = data.frame(signal = s))))
}

#' Quadratic fit: y = a + b*x + c*x^2
fit_quadratic <- function(signal, conc) {
  stopifnot(is.numeric(signal), is.numeric(conc), length(signal) == length(conc))
  model <- lm(conc ~ poly(signal, 2, raw = TRUE))
  list(model = model, name = "quadratic",
       predict = function(s) as.numeric(predict(model, newdata = data.frame(signal = s))))
}

#' Cubic fit: y = a + b*x + c*x^2 + d*x^3
fit_cubic <- function(signal, conc) {
  stopifnot(is.numeric(signal), is.numeric(conc), length(signal) == length(conc))
  model <- lm(conc ~ poly(signal, 3, raw = TRUE))
  list(model = model, name = "cubic",
       predict = function(s) as.numeric(predict(model, newdata = data.frame(signal = s))))
}

#' Logarithmic fit: y = a + b*ln(x)
fit_log <- function(signal, conc) {
  stopifnot(is.numeric(signal), is.numeric(conc), length(signal) == length(conc))
  if (any(signal <= 0)) stop("Log fit requires all signal values > 0")
  model <- lm(conc ~ log(signal))
  list(model = model, name = "log",
       predict = function(s) {
         if (any(s <= 0)) stop("Log fit requires prediction signal values > 0")
         as.numeric(predict(model, newdata = data.frame(signal = s)))
       })
}

#' Power fit: y = a * x^b (linearized via log-log)
fit_power <- function(signal, conc) {
  stopifnot(is.numeric(signal), is.numeric(conc), length(signal) == length(conc))
  if (any(signal <= 0) || any(conc <= 0))
    stop("Power fit requires all signal and concentration values > 0")
  model <- lm(log(conc) ~ log(signal))
  list(model = model, name = "power",
       predict = function(s) {
         if (any(s <= 0)) stop("Power fit requires prediction signal values > 0")
         exp(as.numeric(predict(model, newdata = data.frame(signal = s))))
       })
}

#' Exponential fit: y = a * exp(b*x) (linearized via log)
fit_exponential <- function(signal, conc) {
  stopifnot(is.numeric(signal), is.numeric(conc), length(signal) == length(conc))
  if (any(conc <= 0)) stop("Exponential fit requires all concentration values > 0")
  model <- lm(log(conc) ~ signal)
  list(model = model, name = "exponential",
       predict = function(s) {
         exp(as.numeric(predict(model, newdata = data.frame(signal = s))))
       })
}

#' Estimate 4PL initial parameters from data
.estimate_4pl_start <- function(signal, conc) {
  ord <- order(signal)
  sig_s <- signal[ord]
  conc_s <- conc[ord]

  increasing <- conc_s[length(conc_s)] > conc_s[1]
  if (increasing) {
    a_init <- conc_s[1]
    d_init <- conc_s[length(conc_s)]
    b_init <- -1.0
  } else {
    a_init <- conc_s[length(conc_s)]
    d_init <- conc_s[1]
    b_init <- 1.0
  }

  mid_conc <- (a_init + d_init) / 2
  idx <- which.min(abs(conc_s - mid_conc))
  c_init <- sig_s[idx]
  if (c_init <= 0) c_init <- stats::median(sig_s)

  list(a = a_init, b = b_init, c = c_init, d = d_init)
}

#' Four-parameter logistic (4PL): y = d + (a - d) / (1 + (x/c)^b)
fit_4pl <- function(signal, conc) {
  stopifnot(is.numeric(signal), is.numeric(conc), length(signal) == length(conc))
  if (length(signal) < 5) stop("4PL requires at least 5 data points")

  st <- .estimate_4pl_start(signal, conc)
  df <- data.frame(signal = signal, conc = conc)

  model <- NULL
  for (b_mult in c(1, -1, 0.5, 2)) {
    st_try <- st
    st_try$b <- st$b * b_mult
    model <- tryCatch(
      suppressWarnings(
        nls(conc ~ d + (a - d) / (1 + (signal / c)^b),
            data = df, start = st_try,
            control = nls.control(maxiter = 500, warnOnly = TRUE))),
      error = function(e) NULL
    )
    if (!is.null(model) && !is.null(model$convInfo) && model$convInfo$isConv) break
    model <- NULL
  }
  if (is.null(model)) stop("4PL fitting failed to converge")

  list(model = model, name = "4pl",
       predict = function(s) as.numeric(predict(model, newdata = data.frame(signal = s))))
}

#' Five-parameter logistic (5PL): y = d + (a - d) / (1 + (x/c)^b)^f
fit_5pl <- function(signal, conc) {
  stopifnot(is.numeric(signal), is.numeric(conc), length(signal) == length(conc))
  if (length(signal) < 6) stop("5PL requires at least 6 data points")

  st <- .estimate_4pl_start(signal, conc)
  df <- data.frame(signal = signal, conc = conc)

  model <- NULL
  for (b_mult in c(1, -1, 0.5, 2)) {
    for (f_init in c(1, 0.5, 2)) {
      st_try <- st
      st_try$b <- st$b * b_mult
      st_try$f <- f_init
      model <- tryCatch(
        suppressWarnings(
          nls(conc ~ d + (a - d) / (1 + (signal / c)^b)^f,
              data = df, start = st_try,
              control = nls.control(maxiter = 500, warnOnly = TRUE))),
        error = function(e) NULL
      )
      if (!is.null(model) && !is.null(model$convInfo) && model$convInfo$isConv) break
      model <- NULL
    }
    if (!is.null(model)) break
  }
  if (is.null(model)) stop("5PL fitting failed to converge")

  list(model = model, name = "5pl",
       predict = function(s) as.numeric(predict(model, newdata = data.frame(signal = s))))
}