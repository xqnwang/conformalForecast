test_that("lagmatrix() shifts each column by its own lag or lead", {
  x <- matrix(1:12, nrow = 4, ncol = 3)

  expect_equal(lagmatrix(x, c(0, 0, 0)), x)
  expect_equal(
    lagmatrix(x, c(0, 1, 2)),
    matrix(c(1, 2, 3, 4, NA, 5, 6, 7, NA, NA, 9, 10), nrow = 4)
  )
  expect_equal(
    lagmatrix(x, c(0, -1, -2)),
    matrix(c(1, 2, 3, 4, 6, 7, 8, NA, 11, 12, NA, NA), nrow = 4)
  )
})

test_that("lagmatrix() blanks a column shifted past its own length", {
  x <- matrix(1:8, nrow = 4, ncol = 2)

  expect_true(all(is.na(lagmatrix(x, c(4, 0))[, 1])))
  expect_true(all(is.na(lagmatrix(x, c(0, -9))[, 2])))
  expect_equal(lagmatrix(x, c(4, 0))[, 2], x[, 2])
})

test_that("lagmatrix() keeps the class and the size of its input", {
  x <- ts(matrix(1:12, nrow = 4, ncol = 3), start = 2020, frequency = 4)
  out <- lagmatrix(x, c(0, 1, 2))

  expect_equal(dim(out), dim(x))
  expect_s3_class(out, class(x), exact = TRUE)
})

test_that("lagmatrix() rejects malformed input", {
  x <- matrix(1:12, nrow = 4, ncol = 3)

  expect_error(lagmatrix(1:4, 1), "matrix")
  expect_error(lagmatrix(x, c(1, 2)), "one element per column")
  expect_error(lagmatrix(x, 1:4), "one element per column")
})
