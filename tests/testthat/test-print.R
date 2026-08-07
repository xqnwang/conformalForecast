test_that("print.cvforecast() shows the call, the fit count and the forward step", {
  fc <- test_cv(n = 60, h = 3, initial = 10)
  no_fwd <- test_cv(n = 60, h = 3, initial = 10, forward = FALSE)

  out <- capture.output(print(fc))
  txt <- paste(out, collapse = "\n")

  expect_match(txt, "Cross-validation")
  expect_match(txt, "Call:")
  expect_match(txt, "fit_times = 51")
  expect_match(txt, "Point Forecast")
  # The header must appear once, not once per component of the object.
  expect_equal(sum(grepl("^Cross-validation", out)), 1)

  expect_false(grepl(
    "Point Forecast",
    paste(capture.output(print(no_fwd)), collapse = "\n")
  ))
})

test_that("print() names the conformal method and its cp_times", {
  fc <- test_cv(n = 100, h = 3, initial = 10)

  labels <- c(scp = "SCP", acp = "ACP", pid = "PID", acmcp = "ACMCP")
  fits <- list(
    scp = scp(fc, ncal = 40, rolling = TRUE),
    acp = acp(fc, ncal = 40, rolling = TRUE),
    pid = pid(fc, ncal = 40, rolling = TRUE, scorecastfun = scorecast_fun),
    acmcp = acmcp(fc, ncal = 40, rolling = TRUE)
  )

  for (method in names(fits)) {
    out <- capture.output(print(fits[[method]]))
    txt <- paste(out, collapse = "\n")

    expect_match(txt, labels[[method]], fixed = TRUE)
    expect_match(txt, "cp_times")
    expect_match(txt, "\\(h=1\\)")
    # No stray cross-validation block, and no empty fit_times line.
    expect_false(grepl("fit_times", txt))
    expect_equal(sum(grepl("Call:", out)), 1)
  }
})

test_that("summary() adds the error measures to what print() shows", {
  fc <- test_cv(n = 100, h = 3, initial = 10)
  cp <- scp(fc, symmetric = TRUE, ncal = 40, rolling = TRUE)

  printed <- paste(capture.output(print(cp)), collapse = "\n")
  summarised <- paste(capture.output(print(summary(cp))), collapse = "\n")

  expect_match(summarised, "error measures")
  expect_match(summarised, "RMSE")
  expect_false(grepl("error measures", printed))

  expect_s3_class(summary(fc), "summary.cvforecast")
  expect_s3_class(summary(cp), "summary.cpforecast")
})

test_that("printing never fails on any supported combination", {
  fc <- test_cv(n = 100, h = 2, level = c(80, 95), initial = 10)

  expect_output(print(coverage(fc, level = 95)), "h=1")
  expect_output(print(width(fc, level = 95)), "h=1")

  for (cp in list(
    scp(fc, ncal = 40, rolling = TRUE),
    scp(fc, ncal = 40, symmetric = TRUE),
    acp(fc, ncal = 40, rolling = TRUE),
    pid(fc, ncal = 40, symmetric = TRUE),
    acmcp(fc, ncal = 40)
  )) {
    expect_output(print(cp))
    expect_output(print(summary(cp)))
  }
})
