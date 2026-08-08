# Split a series so that the same data can be reached either in one go or by
# updating a shorter fit.
split_series <- function(n = 100, new = 5) {
  y <- as.numeric(test_series(n))
  list(
    old = ts(y[seq_len(n - new)]),
    new = y[seq.int(n - new + 1L, n)],
    full = ts(y)
  )
}

test_that("update() extends the series and the forecast matrices", {
  parts <- split_series(n = 100, new = 5)
  fc <- cvforecast(
    parts$old,
    forecastfun = naive_fun,
    h = 3,
    level = 95,
    initial = 10
  )
  cp <- scp(fc, ncal = 40, rolling = TRUE)

  updated <- update(cp, new_data = parts$new, forecastfun = naive_fun)

  expect_s3_class(updated, class(cp), exact = TRUE)
  expect_length(updated$x, length(cp$x) + 5)
  expect_equal(nrow(updated$MEAN), nrow(cp$MEAN) + 5)
  expect_aligned_errors(updated)
  expect_valid_intervals(updated)
})

test_that("update() leaves the intervals it issued earlier untouched", {
  parts <- split_series(n = 100, new = 5)
  fc <- cvforecast(
    parts$old,
    forecastfun = naive_fun,
    h = 3,
    level = 95,
    initial = 10
  )

  fits <- list(
    scp = scp(fc, ncal = 40, rolling = TRUE),
    acp = acp(fc, ncal = 40, rolling = TRUE),
    pid = pid(fc, ncal = 40, rolling = TRUE, scorecastfun = scorecast_fun),
    acmcp = acmcp(fc, ncal = 40, rolling = TRUE)
  )

  for (method in names(fits)) {
    cp <- fits[[method]]
    updated <- update(cp, new_data = parts$new, forecastfun = naive_fun)

    before <- as_mat(cp$LOWER[["95%"]])
    after <- as_mat(updated$LOWER[["95%"]])[
      seq_len(nrow(cp$MEAN)),
      ,
      drop = FALSE
    ]
    keep <- !is.na(before)

    expect_equal(after[keep], before[keep], info = method)
  }
})

test_that("update() reproduces a fit computed from the whole series", {
  # The strongest check in the suite: a warm start has to be indistinguishable
  # from replaying the whole history, for any number of new observations.
  for (new in c(1, 5)) {
    parts <- split_series(n = 100, new = new)
    short <- cvforecast(
      parts$old,
      forecastfun = naive_fun,
      h = 3,
      level = 95,
      initial = 10
    )
    long <- cvforecast(
      parts$full,
      forecastfun = naive_fun,
      h = 3,
      level = 95,
      initial = 10
    )

    # Tuning constants that would otherwise be derived from the data are
    # pinned down, so that the two paths compare the same estimator.
    tuning <- list(
      ncal = 40,
      rolling = TRUE,
      Tg = 200,
      delta = 0.01,
      Csat = 0.4,
      KI = 2
    )
    recipes <- list(
      scp = function(x) scp(x, ncal = 40, rolling = TRUE),
      `scp symmetric` = function(x) {
        scp(x, ncal = 40, rolling = TRUE, symmetric = TRUE)
      },
      acp = function(x) acp(x, ncal = 40, rolling = TRUE),
      `acp symmetric` = function(x) {
        acp(x, ncal = 40, rolling = TRUE, symmetric = TRUE)
      },
      pid = function(x) {
        do.call(pid, c(list(x), tuning, list(scorecast = FALSE)))
      },
      `pid scorecaster` = function(x) {
        do.call(pid, c(list(x), tuning, list(scorecastfun = scorecast_fun)))
      },
      acmcp = function(x) do.call(acmcp, c(list(x), tuning))
    )

    for (method in names(recipes)) {
      build <- recipes[[method]]
      label <- paste(method, "new =", new)
      updated <- update(
        build(short),
        new_data = parts$new,
        forecastfun = naive_fun
      )
      scratch <- build(long)

      expect_equal(as_mat(updated$MEAN), as_mat(scratch$MEAN), info = label)
      expect_equal(updated$cp_times, scratch$cp_times, info = label)
      expect_equal(
        as_mat(updated$LOWER[["95%"]]),
        as_mat(scratch$LOWER[["95%"]]),
        info = label
      )
      expect_equal(
        as_mat(updated$UPPER[["95%"]]),
        as_mat(scratch$UPPER[["95%"]]),
        info = label
      )
      # Every new period gets an interval at every horizon.
      expect_equal(
        colSums(!is.na(as_mat(updated$LOWER[["95%"]]))),
        updated$cp_times,
        info = label
      )
    }
  }
})

test_that("update() applies repeatedly", {
  y <- as.numeric(test_series(100))
  fc <- cvforecast(
    ts(y[1:90]),
    forecastfun = naive_fun,
    h = 3,
    level = 95,
    initial = 10
  )

  for (cp in list(
    scp(fc, ncal = 40, rolling = TRUE),
    acp(fc, ncal = 40, rolling = TRUE)
  )) {
    once <- update(cp, new_data = y[91:100], forecastfun = naive_fun)
    twice <- update(cp, new_data = y[91:95], forecastfun = naive_fun)
    twice <- update(twice, new_data = y[96:100], forecastfun = naive_fun)

    expect_equal(as_mat(twice$MEAN), as_mat(once$MEAN))
    expect_equal(as_mat(twice$LOWER[["95%"]]), as_mat(once$LOWER[["95%"]]))
  }
})

test_that("update() does not depend on the environment of the original call", {
  # The stored call refers to `ncal` and `symmetric` by name; if update()
  # re-evaluated it, the object could only be updated from the frame that
  # built it.
  build <- function(y) {
    local_ncal <- 40
    local_symmetric <- TRUE
    fc <- cvforecast(
      y,
      forecastfun = naive_fun,
      h = 3,
      level = 95,
      initial = 10
    )
    scp(fc, ncal = local_ncal, rolling = TRUE, symmetric = local_symmetric)
  }

  parts <- split_series(n = 100, new = 5)
  cp <- build(parts$old)

  expect_no_error(
    updated <- update(cp, new_data = parts$new, forecastfun = naive_fun)
  )
  expect_true(updated$model$args$symmetric)
  expect_equal(updated$model$args$ncal, 40)
})

test_that("update() keeps the tuning constants frozen at their first value", {
  parts <- split_series(n = 100, new = 5)
  fc <- cvforecast(
    parts$old,
    forecastfun = naive_fun,
    h = 3,
    level = 95,
    initial = 10
  )
  cp <- pid(fc, ncal = 40, rolling = TRUE, scorecast = FALSE)

  updated <- update(cp, new_data = parts$new, forecastfun = naive_fun)

  expect_equal(updated$model$args$Tg, cp$model$args$Tg)
  expect_equal(updated$model$args$KI, cp$model$args$KI)
  expect_equal(updated$model$args$Csat, cp$model$args$Csat)
})

test_that("update() carries new exogenous predictors", {
  y <- as.numeric(test_series(100))
  xreg <- matrix(seq_len(103), ncol = 1)

  fc <- cvforecast(
    ts(y[1:95]),
    forecastfun = xreg_fun,
    h = 3,
    level = 95,
    initial = 20,
    xreg = xreg[1:98, , drop = FALSE]
  )
  cp <- scp(fc, ncal = 40, rolling = TRUE)

  updated <- update(
    cp,
    new_data = y[96:100],
    forecastfun = xreg_fun,
    new_xreg = xreg[99:103, , drop = FALSE]
  )

  expect_length(updated$x, 100)
  expect_aligned_errors(updated)
  expect_valid_intervals(updated)
})

test_that("update() rejects mismatched new data", {
  parts <- split_series(n = 100, new = 5)
  fc <- cvforecast(
    parts$old,
    forecastfun = naive_fun,
    h = 3,
    level = 95,
    initial = 10
  )
  cp <- scp(fc, ncal = 40, rolling = TRUE)

  expect_error(update(cp, forecastfun = naive_fun))
  expect_error(update(cp, new_data = numeric(0), forecastfun = naive_fun))
})
