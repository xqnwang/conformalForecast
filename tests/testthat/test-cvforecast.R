test_that("cvforecast() returns the documented structure", {
  h <- 4
  fc <- test_cv(n = 60, h = h, level = 95, initial = 10)

  expect_s3_class(fc, c("cvforecast", "forecast"), exact = TRUE)
  expect_identical(fc$method, "cvforecast")
  expect_contains(
    names(fc),
    c(
      "x",
      "series",
      "MEAN",
      "ERROR",
      "LOWER",
      "UPPER",
      "level",
      "call",
      "forward",
      "fit_times"
    )
  )
  expect_named(fc$LOWER, "95%")
  expect_named(fc$UPPER, "95%")

  # One column per horizon; with `forward = TRUE` the last origin is y_T, so
  # the matrix runs h periods past the end of the series while ERROR stops
  # where the observations do.
  expect_equal(ncol(fc$MEAN), h)
  expect_equal(end(fc$MEAN)[1], end(fc$x)[1] + h)
  expect_equal(end(fc$ERROR)[1], end(fc$x)[1])
  expect_equal(dim(fc$LOWER[["95%"]]), dim(fc$MEAN))
})

test_that("cvforecast() aligns ERROR with MEAN by target period", {
  expect_aligned_errors(test_cv(n = 60, h = 3, initial = 10))
  expect_aligned_errors(test_cv(n = 60, h = 1, initial = 10))
  expect_aligned_errors(test_cv(n = 60, h = 5, initial = 20, window = 30))

  # The same has to hold on a seasonal time base, where the time index is
  # fractional and matching by position is not enough.
  seasonal <- cvforecast(
    test_series(60, frequency = 12),
    forecastfun = naive_fun,
    h = 3,
    level = 95,
    initial = 24
  )
  expect_aligned_errors(seasonal)
  expect_equal(frequency(seasonal$MEAN), 12)
})

test_that("cvforecast() nests its prediction intervals", {
  fc <- test_cv(n = 60, h = 3, level = c(50, 80, 95), initial = 10)
  point <- as_mat(fc$MEAN)

  expect_named(fc$LOWER, c("50%", "80%", "95%"))
  for (lv in names(fc$LOWER)) {
    expect_true(all(as_mat(fc$LOWER[[lv]]) <= point, na.rm = TRUE))
    expect_true(all(as_mat(fc$UPPER[[lv]]) >= point, na.rm = TRUE))
  }
  expect_true(all(
    as_mat(fc$LOWER[["95%"]]) <= as_mat(fc$LOWER[["50%"]]),
    na.rm = TRUE
  ))
  expect_true(all(
    as_mat(fc$UPPER[["95%"]]) >= as_mat(fc$UPPER[["50%"]]),
    na.rm = TRUE
  ))
})

test_that("cvforecast() normalises and sorts the confidence level", {
  as_pct <- test_cv(n = 40, h = 2, level = c(80, 95), initial = 10)
  as_frac <- test_cv(n = 40, h = 2, level = c(0.8, 0.95), initial = 10)
  unsorted <- test_cv(n = 40, h = 2, level = c(95, 50, 80), initial = 10)

  expect_equal(as_pct$level, c(80, 95))
  expect_equal(as_frac$level, c(80, 95))
  expect_equal(as_frac$LOWER, as_pct$LOWER)
  expect_equal(unsorted$level, c(50, 80, 95))
  expect_named(unsorted$LOWER, c("50%", "80%", "95%"))
})

test_that("cvforecast() demands a usable confidence level", {
  y <- test_series(40)

  expect_error(
    cvforecast(y, forecastfun = naive_fun, h = 2, level = NULL, initial = 10),
    "must not be NULL"
  )
  expect_error(
    cvforecast(y, forecastfun = naive_fun, h = 2, level = 150, initial = 10),
    "out of range"
  )
  expect_error(
    cvforecast(y, forecastfun = naive_fun, h = 2, level = -5, initial = 10),
    "out of range"
  )
})

test_that("cvforecast() only produces a final forecast when forward is TRUE", {
  n <- 60
  initial <- 10
  fwd <- test_cv(n = n, h = 3, initial = initial)
  no_fwd <- test_cv(n = n, h = 3, initial = initial, forward = FALSE)

  expect_true(fwd$forward)
  expect_equal(fwd$fit_times, n - initial + 1)
  expect_contains(names(fwd), c("mean", "lower", "upper", "model"))
  expect_length(fwd$mean, 3)

  expect_false(no_fwd$forward)
  expect_equal(no_fwd$fit_times, n - initial)
  expect_false(any(c("mean", "lower", "upper") %in% names(no_fwd)))
  # One origin fewer, so one row fewer.
  expect_equal(nrow(no_fwd$MEAN), nrow(fwd$MEAN) - 1)
})

test_that("cvforecast() honours a rolling window", {
  # `meanf()` averages its whole training window, so a rolling window and an
  # expanding one cannot agree. The first origin is `max(initial, window)`.
  y <- test_series(60)
  rolling <- cvforecast(
    y,
    forecastfun = mean_fun,
    h = 3,
    level = 95,
    initial = 20,
    window = 20
  )
  expanding <- cvforecast(
    y,
    forecastfun = mean_fun,
    h = 3,
    level = 95,
    initial = 20
  )

  expect_equal(dim(rolling$MEAN), dim(expanding$MEAN))
  expect_false(isTRUE(all.equal(as_mat(rolling$MEAN), as_mat(expanding$MEAN))))
  expect_aligned_errors(rolling)

  narrow <- test_cv(n = 60, h = 3, initial = 30, window = 10)
  expect_equal(start(narrow$MEAN)[1], 31)
  expect_equal(narrow$fit_times, 60 - 30 + 1)
})

test_that("cvforecast() warns only about arguments forecastfun cannot take", {
  y <- test_series(40)
  dots_fun <- function(x, h, level, ...) {
    forecast::naive(x, h = h, level = level)
  }

  expect_warning(
    cvforecast(
      y,
      forecastfun = naive_fun,
      h = 2,
      level = 95,
      initial = 10,
      nonexistent = 1
    ),
    "Unused argument"
  )
  expect_no_warning(
    cvforecast(
      y,
      forecastfun = widen_fun,
      h = 2,
      level = 95,
      initial = 10,
      widen = 1
    )
  )
  # A forecaster that absorbs dots can be given anything.
  expect_no_warning(
    cvforecast(
      y,
      forecastfun = dots_fun,
      h = 2,
      level = 95,
      initial = 10,
      anything = 1
    )
  )
})

test_that("cvforecast() passes dots through to forecastfun", {
  y <- test_series(40)

  narrow <- cvforecast(
    y,
    forecastfun = widen_fun,
    h = 2,
    level = 95,
    initial = 10,
    widen = 1
  )
  wide <- cvforecast(
    y,
    forecastfun = widen_fun,
    h = 2,
    level = 95,
    initial = 10,
    widen = 3
  )
  # Filtering the unknown arguments must not take the known ones with it.
  noisy <- suppressWarnings(
    cvforecast(
      y,
      forecastfun = widen_fun,
      h = 2,
      level = 95,
      initial = 10,
      widen = 3,
      nonexistent = 1
    )
  )

  # `widen` only touches the interval, so the point forecasts must agree and
  # the widths must be in a ratio of exactly three.
  expect_equal(as_mat(wide$MEAN), as_mat(narrow$MEAN))
  narrow_width <- as_mat(narrow$UPPER[["95%"]]) - as_mat(narrow$LOWER[["95%"]])
  wide_width <- as_mat(wide$UPPER[["95%"]]) - as_mat(wide$LOWER[["95%"]])
  expect_true(any(!is.na(narrow_width)))
  expect_equal(wide_width, narrow_width * 3)
  expect_equal(as_mat(noisy$UPPER[["95%"]]), as_mat(wide$UPPER[["95%"]]))
})

test_that("cvforecast() carries exogenous predictors through", {
  y <- test_series(60)
  # `forward = TRUE` needs one row of `xreg` per period plus h more.
  xreg <- matrix(seq_len(63), ncol = 1)

  fc <- cvforecast(
    y,
    forecastfun = xreg_fun,
    h = 3,
    level = 95,
    initial = 20,
    xreg = xreg
  )
  plain <- test_cv(n = 60, h = 3, initial = 20)

  expect_contains(names(fc), "xreg")
  expect_aligned_errors(fc)
  # `xreg_fun()` shifts the forecast for period t by `xreg[t]`, so the offset
  # against a plain naive fit says which rows reached `newxreg`.
  observed <- (as_mat(fc$MEAN) - as_mat(plain$MEAN))
  observed <- observed[!is.na(observed)]
  expect_gt(length(observed), 0)
  expect_equal(observed, round(observed))
  expect_true(all(observed > 0))
})
