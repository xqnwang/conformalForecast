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

test_that("all conformal methods calibrate and update point-only forecasts", {
  y <- test_series(85)
  for (forward in c(TRUE, FALSE)) {
    short <- cvforecast(head(y, 80), point_fun, h = 3, level = NULL,
                        initial = 10, forward = forward)
    full <- cvforecast(y, point_fun, h = 3, level = NULL,
                       initial = 10, forward = forward)
    base <- cvforecast(y, naive_fun, h = 3, level = c(80, 90),
                       initial = 10, forward = forward)
    for (method in c("scp", "acp", "pid", "acmcp")) {
      expect_error(conformal(full, method = method), "supply `alpha` explicitly")
      for (symmetric in if (method == "acmcp") FALSE else c(TRUE, FALSE)) {
        args <- list(method = method, alpha = c(0.2, 0.1),
                     ncal = 20, rolling = TRUE)
        if (method != "acmcp") args$symmetric <- symmetric
        if (method %in% c("pid", "acmcp")) {
          args$Csat <- 0.4
          args$KI <- 2
        }
        if (method == "pid") args$scorecastfun <- scorecast_fun
        cp <- do.call(conformal, c(list(object = full), args))
        reference <- do.call(conformal, c(list(object = base), args))
        expect_equal(cp$level, c(80, 90))
        expect_equal(cp$LOWER, reference$LOWER)
        expect_equal(cp$UPPER, reference$UPPER)
        expect_true(any(is.finite(cp$LOWER[["90%"]])))
        expect_valid_intervals(cp, level = "90%")
        expect_equal(coverage(cp, level = 90), coverage(reference, level = 90))
        expect_equal(width(cp, level = 90), width(reference, level = 90))
        expect_equal(accuracy(cp), accuracy(reference))
        expect_output(print(summary(cp)), "error measures")
        before <- do.call(conformal, c(list(object = short), args))
        updated <- update(before, new_data = as.numeric(tail(y, 5)),
                          forecastfun = point_fun)
        expect_equal(updated$MEAN, cp$MEAN)
        expect_equal(updated$ERROR, cp$ERROR)
        expect_equal(updated$LOWER, cp$LOWER)
        expect_equal(updated$UPPER, cp$UPPER)
        expect_equal(updated$mean, cp$mean)
        expect_equal(updated$lower, cp$lower)
        expect_equal(updated$upper, cp$upper)
        expect_null(updated$model$cvforecast$args$level)
      }
    }
  }
})

test_that("point-only updates preserve seasonal time and external regressors", {
  y <- test_series(85, frequency = 12)
  f <- function(x, h, level, xreg, newxreg) {
    stopifnot(is.null(level), NROW(xreg) == length(x))
    list(mean = rep(as.numeric(tail(x, 1)), h) + rowSums(newxreg))
  }
  for (forward in c(TRUE, FALSE)) {
    for (h in c(1, 3)) {
      xr <- cbind(a = seq_len(85 + h - !forward) / 100, b = 0.1)
      build <- function(n) {
        cvforecast(window(y, end = time(y)[n]), f, h = h, level = NULL,
                   window = 20, forward = forward,
                   xreg = xr[seq_len(n + h - !forward), , drop = FALSE])
      }
      short <- build(80)
      full <- build(85)
      for (method in c("scp", "acp", "pid", "acmcp")) {
        args <- list(method = method, alpha = 0.1, ncal = 20, rolling = TRUE)
        if (method %in% c("pid", "acmcp")) {
          args$Csat <- 0.4
          args$KI <- 2
        }
        if (method == "pid") args$scorecastfun <- scorecast_fun
        updated <- do.call(conformal, c(list(object = short), args))
        for (rows in list(81:82, 83:85)) {
          updated <- update(updated, new_data = as.numeric(y[rows]),
                            forecastfun = f,
                            new_xreg = xr[rows + h - !forward, , drop = FALSE])
        }
        expected <- do.call(conformal, c(list(object = full), args))
        expect_equal(updated$MEAN, expected$MEAN)
        expect_equal(updated$ERROR, expected$ERROR)
        expect_equal(updated$LOWER, expected$LOWER)
        expect_equal(updated$UPPER, expected$UPPER)
        expect_equal(updated$mean, expected$mean)
        expect_equal(updated$lower, expected$lower)
        expect_equal(updated$upper, expected$upper)
        expect_equal(updated$xreg, expected$xreg)
      }
    }
  }
})
