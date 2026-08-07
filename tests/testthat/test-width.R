test_that("width() is the distance between the two bounds", {
  fc <- test_cv(n = 100, h = 3, initial = 10)
  wid <- width(fc, level = 95)

  expect_s3_class(wid, "width")
  expect_named(wid, c("width", "mean"))
  expect_named(wid$mean, c("h=1", "h=2", "h=3"))

  expected <- as_mat(fc$UPPER[["95%"]]) - as_mat(fc$LOWER[["95%"]])
  expect_equal(as_mat(wid$width), expected)
  expect_true(all(as_mat(wid$width) >= 0, na.rm = TRUE))
  expect_equal(
    as.numeric(wid$mean),
    as.numeric(colMeans(expected, na.rm = TRUE))
  )
})

test_that("width() adds the median and the rolling summaries only when asked", {
  fc <- test_cv(n = 100, h = 3, initial = 10)

  plain <- width(fc, level = 95)
  rolled <- width(fc, level = 95, includemedian = TRUE, window = 20)

  expect_false(any(c("median", "rollmean") %in% names(plain)))
  expect_contains(names(rolled), c("median", "rollmean", "rollmedian"))
  expect_equal(
    as.numeric(rolled$median),
    as.numeric(apply(as_mat(rolled$width), 2, median, na.rm = TRUE))
  )
  expect_equal(ncol(rolled$rollmean), 3)
  expect_equal(dim(rolled$rollmedian), dim(rolled$rollmean))
  expect_true(all(rolled$rollmean >= 0, na.rm = TRUE))
})

test_that("width() accepts the pieces instead of the object", {
  fc <- test_cv(n = 100, h = 3, initial = 10)

  from_object <- width(fc, level = 95, includemedian = TRUE, window = 20)
  from_parts <- width(
    LOWER = fc$LOWER,
    UPPER = fc$UPPER,
    level = 95,
    includemedian = TRUE,
    window = 20
  )

  expect_equal(from_parts, from_object)
})

test_that("width() refuses what it cannot compute", {
  fc <- test_cv(n = 60, h = 2, level = 95, initial = 10)

  expect_error(width(fc, level = 80), "confidence level")
  expect_error(width(fc, level = c(80, 95)), "one confidence level")
  expect_error(width(fc, level = 150), "out of range")
  expect_error(width(LOWER = fc$LOWER, level = 95), "required")
})
