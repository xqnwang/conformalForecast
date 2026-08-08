# This file holds the contract every conformal method has to honour, so that
# the per-method files only have to cover what is specific to them.

conformal_fits <- function(fc) {
  list(
    scp = scp(fc, ncal = 40, rolling = TRUE),
    acp = acp(fc, ncal = 40, rolling = TRUE),
    pid = pid(fc, ncal = 40, rolling = TRUE, scorecastfun = scorecast_fun),
    acmcp = acmcp(fc, ncal = 40, rolling = TRUE)
  )
}

test_that("conformal() dispatches to the method it is asked for", {
  fc <- test_cv(n = 100, h = 3, initial = 10)

  expect_equal(
    as_mat(conformal(fc, method = "scp", ncal = 40, rolling = TRUE)$LOWER[[
      "95%"
    ]]),
    as_mat(scp(fc, ncal = 40, rolling = TRUE)$LOWER[["95%"]])
  )
  expect_equal(
    as_mat(conformal(fc, method = "acp", ncal = 40, rolling = TRUE)$LOWER[[
      "95%"
    ]]),
    as_mat(acp(fc, ncal = 40, rolling = TRUE)$LOWER[["95%"]])
  )
  expect_equal(
    as_mat(conformal(fc, method = "pid", ncal = 40, scorecast = FALSE)$LOWER[[
      "95%"
    ]]),
    as_mat(pid(fc, ncal = 40, scorecast = FALSE)$LOWER[["95%"]])
  )
  expect_equal(
    as_mat(conformal(fc, method = "acmcp", ncal = 40)$LOWER[["95%"]]),
    as_mat(acmcp(fc, ncal = 40)$LOWER[["95%"]])
  )
  # The first option of `method` is the default.
  expect_s3_class(conformal(fc, ncal = 40, rolling = TRUE), "scp")
})

test_that("conformal() rejects what it cannot dispatch on", {
  fc <- test_cv(n = 60, h = 2, initial = 10)

  expect_error(conformal(fc, method = "nonesuch", ncal = 20))
  expect_error(conformal(forecast::naive(test_series(30), h = 2)))
})

test_that("every conformal method honours the shared contract", {
  fc <- test_cv(n = 100, h = 3, level = c(80, 95), initial = 10)
  fits <- conformal_fits(fc)

  for (method in names(fits)) {
    cp <- fits[[method]]

    expect_cpforecast(cp, method)
    expect_valid_intervals(cp, "80%")
    # The point forecasts and errors belong to the cross-validation, not to
    # the conformal step, so no method may touch them.
    expect_equal(as_mat(cp$MEAN), as_mat(fc$MEAN), info = method)
    expect_equal(as_mat(cp$ERROR), as_mat(fc$ERROR), info = method)
    expect_equal(cp$level, fc$level, info = method)
    expect_equal(cp$x, fc$x, info = method)
    # One interval per horizon per origin from `ncal + h - 1` onwards.
    expect_equal(
      colSums(!is.na(as_mat(cp$LOWER[["95%"]]))),
      cp$cp_times,
      info = method
    )
    # The forward step is carried over from the cross-validation.
    expect_equal(cp$mean, fc$mean, info = method)
    # The resolved arguments are recorded so that update() can replay them.
    expect_contains(names(cp$model$args), c("alpha", "ncal", "rolling"))
    expect_contains(names(cp$model$cvforecast), c("call", "args"))
  }
})

test_that("every conformal method works without the forward step", {
  fc <- test_cv(n = 100, h = 3, initial = 10, forward = FALSE)
  fits <- conformal_fits(fc)

  for (method in names(fits)) {
    expect_cpforecast(fits[[method]], method)
    expect_false(any(c("mean", "lower", "upper") %in% names(fits[[method]])))
  }
})

test_that("symmetric scores give intervals centred on the point forecast", {
  fc <- test_cv(n = 100, h = 3, initial = 10)

  # `acmcp()` is defined on signed errors and has no `symmetric` argument.
  for (cp in list(
    scp(fc, symmetric = TRUE, ncal = 40, rolling = TRUE),
    acp(fc, symmetric = TRUE, ncal = 40, rolling = TRUE),
    pid(fc, symmetric = TRUE, ncal = 40, rolling = TRUE)
  )) {
    point <- as_mat(cp$MEAN)
    expect_equal(
      point - as_mat(cp$LOWER[["95%"]]),
      as_mat(cp$UPPER[["95%"]]) - point
    )
  }
})

test_that("every conformal method reaches roughly the nominal coverage", {
  fc <- test_cv(n = 150, h = 3, level = 80, initial = 10)
  fits <- conformal_fits(fc)

  for (method in names(fits)) {
    observed <- coverage(fits[[method]], level = 80)$mean
    expect_true(all(observed > 0.6), info = method)
    expect_true(all(observed <= 1), info = method)
  }
})

test_that("scp nests the intervals of nested levels", {
  fc <- test_cv(n = 100, h = 3, level = c(80, 95), initial = 10)
  cp <- scp(fc, ncal = 40, rolling = TRUE)

  expect_true(all(
    as_mat(cp$UPPER[["95%"]]) - as_mat(cp$LOWER[["95%"]]) >=
      as_mat(cp$UPPER[["80%"]]) - as_mat(cp$LOWER[["80%"]]),
    na.rm = TRUE
  ))
})

test_that("every conformal method rejects the settings it cannot work with", {
  fc <- test_cv(n = 60, h = 2, initial = 10)
  methods <- list(
    scp = scp,
    acp = acp,
    pid = function(...) pid(..., scorecast = FALSE),
    acmcp = acmcp
  )

  for (method in names(methods)) {
    f <- methods[[method]]
    expect_error(f(fc, alpha = 0, ncal = 20), "alpha", info = method)
    expect_error(f(fc, alpha = 1.5, ncal = 20), "alpha", info = method)
    expect_error(f(fc, ncal = 1000), "ncal", info = method)
  }
})
