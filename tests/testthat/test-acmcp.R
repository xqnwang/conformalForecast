test_that("acmcp() produces asymmetric intervals only", {
  fc <- test_cv(n = 100, h = 3, initial = 10)
  cp <- acmcp(fc, ncal = 40, rolling = TRUE)

  # There is no `symmetric` argument: the method is defined on signed errors.
  expect_false("symmetric" %in% names(formals(acmcp)))
  expect_false("symmetric" %in% names(cp$model$args))

  point <- as_mat(cp$MEAN)
  expect_false(isTRUE(all.equal(
    point - as_mat(cp$LOWER[["95%"]]),
    as_mat(cp$UPPER[["95%"]]) - point
  )))
})

test_that("acmcp() resolves its tuning constants and stores its state", {
  fc <- test_cv(n = 100, h = 3, initial = 10)
  cp <- acmcp(fc, ncal = 40)

  expect_equal(cp$model$args$Tg, NROW(fc$ERROR))
  expect_equal(cp$model$args$KI, max(abs(fc$ERROR), na.rm = TRUE))
  expect_equal(
    cp$model$args$Csat,
    2 / pi * (ceiling(log(cp$model$args$Tg) * 0.01) - 1 / log(cp$model$args$Tg))
  )
  expect_contains(
    names(cp$model),
    c("args", "lr_update", "t_last", "state", "scorecaster")
  )
  expect_equal(ncol(cp$model$lr_update), ncol(cp$MEAN))
})

test_that("acmcp() modules each change the answer", {
  fc <- test_cv(n = 100, h = 3, initial = 10)
  base <- acmcp(
    fc,
    ncal = 40,
    rolling = TRUE,
    integrate = FALSE,
    scorecast = FALSE
  )

  variants <- list(
    integrator = acmcp(
      fc,
      ncal = 40,
      rolling = TRUE,
      integrate = TRUE,
      scorecast = FALSE
    ),
    scorecaster = acmcp(
      fc,
      ncal = 40,
      rolling = TRUE,
      integrate = FALSE,
      scorecast = TRUE
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
