# Fixtures ---------------------------------------------------------------

# A smooth deterministic series. No RNG is involved, so the tests do not
# depend on the seed, on the RNG version, or on the order in which they run.
test_series <- function(n = 100, frequency = 1L) {
  t <- seq_len(n)
  ts(10 + 0.05 * t + sin(t / 3) + 0.3 * cos(t / 7), frequency = frequency)
}

# `naive()` costs almost nothing to fit, so a full cross-validation over 100
# observations stays well under a tenth of a second.
naive_fun <- function(x, h, level) {
  forecast::naive(x, h = h, level = level)
}

scorecast_fun <- function(x, h) {
  forecast::naive(x, h = h)
}

# Sensitive to the width of the training window, unlike `naive()`, which only
# ever looks at the last observation.
mean_fun <- function(x, h, level) {
  forecast::meanf(x, h = h, level = level)
}

# A forecaster with an extra argument of its own, used to check that `...`
# reaches `forecastfun`.
widen_fun <- function(x, h, level, widen) {
  fit <- forecast::naive(x, h = h, level = level)
  centre <- as.numeric(fit$mean)
  fit$lower[] <- centre - widen * (centre - fit$lower)
  fit$upper[] <- centre + widen * (fit$upper - centre)
  fit
}

# A forecaster that uses an exogenous predictor without fitting anything, so
# that the tests of the xreg plumbing cannot fail for numerical reasons. The
# argument names are the ones cvforecast() documents.
xreg_fun <- function(x, h, level, xreg, newxreg) {
  fit <- forecast::naive(x, h = h, level = level)
  # A one-column predictor arrives as a plain vector.
  shift <- as.numeric(as.matrix(newxreg)[, 1])
  fit$mean <- fit$mean + shift
  fit$lower[] <- fit$lower + shift
  fit$upper[] <- fit$upper + shift
  fit
}

# `ncal = 40` is the smallest calibration set that keeps the asymmetric 95%
# quantile finite: the conformal quantile is `Inf` whenever
# `ceiling((ncal + 1) * (1 - alpha / 2)) > ncal`.
test_cv <- function(n = 100, h = 3, level = 95, initial = 10, ...) {
  cvforecast(
    test_series(n),
    forecastfun = naive_fun,
    h = h,
    level = level,
    initial = initial,
    ...
  )
}

# Drop the ts attributes so that matrices with different time bases can be
# compared cell by cell.
as_mat <- function(x) {
  matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x))
}

# Expectations -----------------------------------------------------------

# The two bounds are not checked for ordering: `pid()` tracks each bound with
# its own recursion, so the tracked quantile can go negative and the interval
# can collapse or invert. That is a property of the method, not a defect.
expect_valid_intervals <- function(object, level = "95%") {
  lower <- as_mat(object$LOWER[[level]])
  upper <- as_mat(object$UPPER[[level]])
  point <- as_mat(object$MEAN)

  expect_equal(dim(lower), dim(point))
  expect_equal(dim(upper), dim(point))
  expect_identical(is.na(lower), is.na(upper))

  invisible(object)
}

expect_cpforecast <- function(object, method) {
  expect_s3_class(
    object,
    c(method, "cpforecast", "cvforecast", "forecast"),
    exact = TRUE
  )
  expect_identical(object$method, method)
  expect_length(object$cp_times, ncol(object$MEAN))
  expect_true(all(diff(object$cp_times) < 0))
  expect_valid_intervals(object)

  invisible(object)
}

# `ERROR[t, h]` must be exactly `x[t] - MEAN[t, h]`. This is the invariant
# that every conformal method depends on, and it is the one that breaks first
# if the horizon alignment is off by one.
expect_aligned_errors <- function(object) {
  errors <- object$ERROR
  point <- stats::window(object$MEAN, start = start(errors), end = end(errors))
  i <- match(
    round(as.numeric(time(errors)), 8),
    round(as.numeric(time(object$x)), 8)
  )
  actual <- as.numeric(object$x)[i]

  expect_false(anyNA(i))
  expect_equal(as_mat(errors), actual - as_mat(point))

  invisible(object)
}
