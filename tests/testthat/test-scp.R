test_that("scp() replaces the intervals it was given", {
  fc <- test_cv(n = 100, h = 3, initial = 10)
  cp <- scp(fc, ncal = 40, rolling = TRUE)

  expect_false(isTRUE(all.equal(
    as_mat(cp$LOWER[["95%"]]),
    as_mat(fc$LOWER[["95%"]])
  )))
  # Asymmetric scores are the default and give a one-sided correction.
  point <- as_mat(cp$MEAN)
  expect_false(isTRUE(all.equal(
    point - as_mat(cp$LOWER[["95%"]]),
    as_mat(cp$UPPER[["95%"]]) - point
  )))
})

test_that("scp() distinguishes rolling from expanding calibration sets", {
  fc <- test_cv(n = 100, h = 3, initial = 10)

  rolling <- scp(fc, ncal = 40, rolling = TRUE)
  expanding <- scp(fc, ncal = 40, rolling = FALSE)

  expect_false(isTRUE(all.equal(
    as_mat(rolling$LOWER[["95%"]]),
    as_mat(expanding$LOWER[["95%"]])
  )))
  # Both start conformalising at the same point.
  expect_equal(rolling$cp_times, expanding$cp_times)
  expect_equal(diff(rolling$cp_times), c(-1, -1))
})

test_that("scp() treats constant weights as no weights", {
  fc <- test_cv(n = 100, h = 3, initial = 10)

  unweighted <- scp(fc, ncal = 40, rolling = TRUE)
  weighted <- scp(fc, ncal = 40, rolling = TRUE, weightfun = function(n) {
    rep(1, n)
  })

  expect_equal(
    as_mat(weighted$LOWER[["95%"]]),
    as_mat(unweighted$LOWER[["95%"]])
  )
})

test_that("scp() reacts to a weight function and to Kish's sample size", {
  # `kess` rescales the sample size, which only moves an interpolating
  # estimator: quantile types 1 to 3 are step functions and do not react.
  fc <- test_cv(n = 120, h = 3, level = 80, initial = 10)
  decay <- function(n) 0.999^(n:1)

  unweighted <- scp(
    fc,
    ncal = 60,
    rolling = TRUE,
    symmetric = TRUE,
    quantiletype = 7
  )
  weighted <- scp(
    fc,
    ncal = 60,
    rolling = TRUE,
    symmetric = TRUE,
    weightfun = decay,
    quantiletype = 7
  )
  effective <- scp(
    fc,
    ncal = 60,
    rolling = TRUE,
    symmetric = TRUE,
    weightfun = decay,
    quantiletype = 7,
    kess = TRUE
  )

  expect_false(isTRUE(all.equal(
    as_mat(weighted$LOWER[["80%"]]),
    as_mat(unweighted$LOWER[["80%"]])
  )))
  expect_false(isTRUE(all.equal(
    as_mat(effective$LOWER[["80%"]]),
    as_mat(weighted$LOWER[["80%"]])
  )))
  expect_valid_intervals(effective, "80%")
})

test_that("scp() records the arguments it actually used", {
  fc <- test_cv(n = 100, h = 3, initial = 10)
  cp <- scp(fc, ncal = 40, rolling = TRUE, symmetric = TRUE)

  expect_identical(cp$model$args$ncal, 40)
  expect_true(cp$model$args$rolling)
  expect_true(cp$model$args$symmetric)
  expect_equal(cp$model$args$quantiletype, 1)
})
