test_that("coverage() returns a proportion per horizon", {
  fc <- test_cv(n = 100, h = 3, initial = 10)
  cov <- coverage(fc, level = 95)

  expect_s3_class(cov, "coverage")
  expect_named(cov, c("mean", "ifinn"))
  expect_named(cov$mean, c("h=1", "h=2", "h=3"))
  expect_equal(ncol(cov$ifinn), 3)
  expect_type(cov$ifinn, "logical")
  expect_true(all(cov$mean >= 0 & cov$mean <= 1))
  # The mean is the mean of the indicator matrix, column by column.
  expect_equal(
    as.numeric(cov$mean),
    as.numeric(colMeans(as_mat(cov$ifinn), na.rm = TRUE))
  )
})

test_that("coverage() flags the periods that fall outside the interval", {
  fc <- test_cv(n = 100, h = 3, initial = 10)
  cov <- coverage(fc, level = 95)

  errors <- fc$ERROR
  point <- stats::window(fc$MEAN, start = start(errors), end = end(errors))
  lower <- stats::window(
    fc$LOWER[["95%"]],
    start = start(errors),
    end = end(errors)
  )
  upper <- stats::window(
    fc$UPPER[["95%"]],
    start = start(errors),
    end = end(errors)
  )
  actual <- as_mat(point) + as_mat(errors)

  inside <- (actual >= as_mat(lower)) & (actual <= as_mat(upper))
  expect_equal(as_mat(cov$ifinn), inside * 1)
})

test_that("coverage() adds a rolling mean when asked", {
  fc <- test_cv(n = 100, h = 3, initial = 10)

  plain <- coverage(fc, level = 95)
  rolled <- coverage(fc, level = 95, window = 20)

  expect_false("rollmean" %in% names(plain))
  expect_named(rolled, c("mean", "ifinn", "rollmean"))
  expect_equal(ncol(rolled$rollmean), 3)
  expect_lt(nrow(rolled$rollmean), nrow(rolled$ifinn))
  expect_true(all(rolled$rollmean >= 0 & rolled$rollmean <= 1, na.rm = TRUE))
})

test_that("coverage() accepts the pieces instead of the object", {
  fc <- test_cv(n = 100, h = 3, initial = 10)

  from_object <- coverage(fc, level = 95, window = 20)
  from_parts <- coverage(
    x = fc$x,
    LOWER = fc$LOWER,
    UPPER = fc$UPPER,
    level = 95,
    window = 20
  )
  # A single matrix per bound is accepted as well as the list of levels.
  from_matrix <- coverage(
    x = fc$x,
    LOWER = fc$LOWER[["95%"]],
    UPPER = fc$UPPER[["95%"]],
    level = 95,
    window = 20
  )

  expect_equal(from_parts, from_object)
  expect_equal(from_matrix, from_object)
})

test_that("coverage() refuses what it cannot compute", {
  fc <- test_cv(n = 60, h = 2, level = 95, initial = 10)

  expect_error(coverage(fc, level = 80), "confidence level")
  expect_error(coverage(fc, level = c(80, 95)), "one confidence level")
  expect_error(coverage(fc, level = 150), "out of range")
  expect_error(coverage(x = fc$x, LOWER = fc$LOWER, level = 95), "required")
})
