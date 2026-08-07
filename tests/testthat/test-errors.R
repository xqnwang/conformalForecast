test_that("the point measures match their definitions", {
  resid <- c(1, -3, 2)
  actual <- c(10, 20, 40)

  expect_equal(ME(resid), 0)
  expect_equal(MAE(resid), 2)
  expect_equal(MSE(resid), 14 / 3)
  expect_equal(RMSE(resid), sqrt(14 / 3))
  expect_equal(MPE(resid, actual), mean(c(10, -15, 5)))
  expect_equal(MAPE(resid, actual), mean(c(10, 15, 5)))

  # RMSE is the root of a mean, so it can never fall below MAE, and the two
  # coincide only when every residual has the same magnitude.
  expect_gte(RMSE(resid), MAE(resid))
  expect_equal(RMSE(c(2, -2, 2)), MAE(c(2, -2, 2)))
})

test_that("the scaled measures divide by the training scale", {
  resid <- c(1, -3, 2)
  train <- c(1, 3, 2, 6, 4)

  expect_equal(
    MASE(resid, train = train, period = 1),
    mean(abs(resid)) / mean(abs(diff(train)))
  )
  expect_equal(
    RMSSE(resid, train = train, period = 1),
    sqrt(mean(resid^2 / mean(diff(train)^2)))
  )

  # Seasonal data is differenced at the seasonal lag instead.
  seasonal_train <- c(1, 3, 2, 6, 4, 8, 5, 10)
  expect_equal(
    MASE(resid, train = seasonal_train, period = 4),
    mean(abs(resid)) / mean(abs(diff(seasonal_train, lag = 4)))
  )
})

test_that("winkler_score() charges the width plus the distance of a miss", {
  alpha <- 0.05

  inside <- winkler_score(
    lower = c(0, 1),
    upper = c(2, 3),
    actual = c(1, 2),
    level = 95
  )
  below <- winkler_score(lower = 0, upper = 2, actual = -1, level = 95)
  above <- winkler_score(lower = 0, upper = 2, actual = 3, level = 95)

  expect_equal(inside, 2)
  expect_equal(below, 2 + (2 / alpha) * 1)
  expect_equal(above, 2 + (2 / alpha) * 1)
  # The level may be given as a fraction.
  expect_equal(
    winkler_score(lower = 0, upper = 2, actual = -1, level = 0.95),
    below
  )
  expect_error(winkler_score(0, 2, 1, level = 150), "out of range")
})

test_that("MSIS() is the Winkler score on the training scale", {
  train <- c(1, 3, 2, 6, 4)

  expect_equal(
    MSIS(
      lower = 0,
      upper = 2,
      actual = 1,
      train = train,
      level = 95,
      period = 1
    ),
    2 / mean(abs(diff(train)))
  )
})

test_that("the measure lists expose the exported functions", {
  expect_named(
    point_measures,
    c("ME", "MAE", "MSE", "RMSE", "MPE", "MAPE", "MASE", "RMSSE")
  )
  expect_named(interval_measures, c("Winkler", "MSIS"))
  expect_true(all(vapply(point_measures, is.function, logical(1))))
  expect_true(all(vapply(interval_measures, is.function, logical(1))))
})

test_that("the measures drop missing values by default", {
  resid <- c(1, NA, -3, 2)

  expect_equal(MAE(resid), 2)
  expect_equal(MAE(resid, na.rm = FALSE), NA_real_)
  expect_equal(ME(resid), 0)
})

test_that("accuracy() returns one row per requested horizon", {
  fc <- test_cv(n = 100, h = 3, initial = 10)
  cp <- scp(fc, symmetric = TRUE, ncal = 40, rolling = TRUE)

  pooled <- accuracy(cp, measures = point_measures)
  by_h <- accuracy(cp, measures = point_measures, byhorizon = TRUE)
  intervals <- accuracy(cp, measures = interval_measures, level = 95)

  expect_true(is.matrix(pooled))
  expect_equal(rownames(pooled), "CV")
  expect_equal(rownames(by_h), c("CV h=1", "CV h=2", "CV h=3"))
  expect_equal(colnames(by_h), names(point_measures))
  expect_equal(colnames(intervals), c("Winkler_95", "MSIS_95"))
  expect_true(all(is.finite(intervals)) && all(intervals > 0))
})

test_that("accuracy() works on a plain cross-validation object", {
  fc <- test_cv(n = 100, h = 3, initial = 10)

  out <- accuracy(fc, measures = point_measures, byhorizon = TRUE)
  with_test_set <- accuracy(fc, x = c(13, 13.5, 14), measures = point_measures)

  expect_equal(nrow(out), 3)
  # Errors accumulate with the horizon for a naive forecaster on a trending
  # series.
  expect_true(all(diff(out[, "MAE"]) > 0))
  expect_gt(nrow(with_test_set), 1)
  expect_contains(rownames(with_test_set), "CV")
})

test_that("accuracy() rejects a period that contradicts the series", {
  fc <- cvforecast(
    test_series(100, frequency = 4),
    forecastfun = naive_fun,
    h = 3,
    level = 95,
    initial = 12
  )

  expect_error(
    accuracy(fc, measures = point_measures, period = 7),
    "does not match"
  )
})
