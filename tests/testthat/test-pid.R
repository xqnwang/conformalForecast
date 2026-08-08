test_that("pid() needs a scorecaster whenever scorecasting is on", {
  fc <- test_cv(n = 60, h = 2, initial = 10)

  expect_error(
    pid(fc, ncal = 20, scorecast = TRUE, scorecastfun = NULL),
    "scorecastfun"
  )
  # The default is `!symmetric`, so an asymmetric call without a scorecaster
  # fails while a symmetric one does not.
  expect_error(pid(fc, ncal = 20, symmetric = FALSE), "scorecastfun")
  expect_no_error(pid(fc, ncal = 20, symmetric = TRUE))
})

test_that("pid() resolves and records its tuning constants", {
  fc <- test_cv(n = 100, h = 3, initial = 10)

  derived <- pid(fc, ncal = 40, scorecast = FALSE)
  expect_equal(derived$model$args$Tg, NROW(fc$ERROR))
  expect_equal(derived$model$args$KI, max(abs(fc$ERROR), na.rm = TRUE))
  expect_equal(derived$model$args$delta, 0.01)
  expect_equal(
    derived$model$args$Csat,
    2 /
      pi *
      (ceiling(log(derived$model$args$Tg) * 0.01) -
        1 / log(derived$model$args$Tg))
  )
  expect_gt(derived$model$args$Csat, 0)

  given <- pid(
    fc,
    ncal = 40,
    scorecast = FALSE,
    Tg = 200,
    delta = 0.02,
    Csat = 0.4,
    KI = 2,
    lr = 0.05
  )
  expect_equal(
    given$model$args[c("Tg", "delta", "Csat", "KI", "lr")],
    list(Tg = 200, delta = 0.02, Csat = 0.4, KI = 2, lr = 0.05)
  )
})

test_that("pid() modules each change the answer", {
  fc <- test_cv(n = 100, h = 3, initial = 10)
  base <- pid(
    fc,
    ncal = 40,
    rolling = TRUE,
    scorecast = FALSE,
    integrate = FALSE
  )

  variants <- list(
    integrator = pid(
      fc,
      ncal = 40,
      rolling = TRUE,
      scorecast = FALSE,
      integrate = TRUE
    ),
    scorecaster = pid(
      fc,
      ncal = 40,
      rolling = TRUE,
      integrate = FALSE,
      scorecastfun = scorecast_fun
    ),
    `learning rate` = pid(
      fc,
      ncal = 40,
      rolling = TRUE,
      scorecast = FALSE,
      integrate = FALSE,
      lr = 1
    )
  )

  for (name in names(variants)) {
    expect_false(
      isTRUE(all.equal(
        as_mat(variants[[name]]$LOWER[["95%"]]),
        as_mat(base$LOWER[["95%"]])
      )),
      info = name
    )
  }
})

test_that("pid() stores the state a warm start needs", {
  fc <- test_cv(n = 100, h = 3, initial = 10)
  cp <- pid(fc, ncal = 40, rolling = TRUE, scorecast = FALSE)

  expect_contains(names(cp$model), c("args", "lr_update", "t_last", "state"))
  expect_equal(ncol(cp$model$lr_update), ncol(cp$MEAN))
  expect_length(cp$model$state, ncol(cp$MEAN))
  # `t_last` marks how far the recursion has already been run.
  expect_length(cp$model$t_last, 1)
  expect_gt(cp$model$t_last, 0)
})
