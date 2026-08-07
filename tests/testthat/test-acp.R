test_that("acp() adapts, and gamma controls how much", {
  fc <- test_cv(n = 100, h = 3, initial = 10)

  slow <- acp(fc, gamma = 0.001, ncal = 40, rolling = TRUE)
  fast <- acp(fc, gamma = 0.05, ncal = 40, rolling = TRUE)

  expect_false(isTRUE(all.equal(
    as_mat(slow$LOWER[["95%"]]),
    as_mat(fast$LOWER[["95%"]])
  )))
  # A larger step size moves alpha_t further from its starting point. The
  # updating path is kept in `model$alpha_update`, one `n x h` matrix of
  # alpha_t values per nominal level, split into `lower`/`upper` when the
  # method is asymmetric.
  slow_range <- diff(range(
    slow$model$alpha_update$lower[["95%"]],
    na.rm = TRUE
  ))
  fast_range <- diff(range(
    fast$model$alpha_update$lower[["95%"]],
    na.rm = TRUE
  ))
  expect_gt(fast_range, slow_range)
})

test_that("acp() stores the alpha updating path", {
  # Unlike pid() and acmcp(), which carry a compact `state`/`t_last` warm
  # start, acp() records the whole alpha_t path and update() replays it.
  # `alpha_update` therefore has to hold, for every nominal level, one
  # `n x h` matrix of the alpha_t actually used at each origin and horizon.
  fc <- test_cv(n = 100, h = 3, initial = 10)
  levels <- names(fc$LOWER)
  shape <- c(nrow(fc$MEAN), ncol(fc$MEAN))

  asym <- acp(fc, gamma = 0.01, ncal = 40, rolling = TRUE)
  expect_contains(names(asym$model), c("args", "alpha_update"))
  expect_named(asym$model$alpha_update, c("lower", "upper"))
  for (bound in c("lower", "upper")) {
    path <- asym$model$alpha_update[[bound]]
    expect_named(path, levels, info = bound)
    expect_equal(dim(path[[1]]), shape, info = bound)
    # The path is only defined once calibration can start, and one origin
    # later for each extra horizon, so it is padded with NAs at the front.
    expect_equal(
      colSums(!is.na(as_mat(path[[1]]))),
      asym$cp_times,
      info = bound
    )
  }
  expect_identical(asym$model$args$gamma, 0.01)

  sym <- acp(fc, ncal = 40, rolling = TRUE, symmetric = TRUE)
  expect_named(sym$model$alpha_update, "alpha")
  expect_named(sym$model$alpha_update$alpha, levels)
  expect_equal(dim(sym$model$alpha_update$alpha[[1]]), shape)
})

test_that("acp() rejects a step size it cannot use", {
  fc <- test_cv(n = 60, h = 2, initial = 10)

  expect_error(acp(fc, gamma = -1, ncal = 20), "gamma")
})
