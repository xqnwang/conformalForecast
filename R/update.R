#' Update and reperform cross-validation forecasting and conformal prediction
#'
#' Update conformal prediction intervals and other information by applying the
#' \code{cvforecast} and \code{conformal} functions.
#'
#' @param object An object of class \code{"cpforecast"}.
#' @param new_data A non-empty numeric vector of newly available data.
#' @param forecastfun Function to return an object of class \code{"forecast"}.
#' Its first argument must be a univariate time series, and it must have an
#' argument \code{h} for the forecast horizon and an argument \code{level} for
#' the confidence level for prediction intervals. If exogenous predictors are used,
#' then it must also have \code{xreg} and \code{newxreg} arguments corresponding
#' to the training and test periods, respectively.
#' @param new_xreg Newly available exogenous predictor variables passed to
#' \code{forecastfun} if required. The number of rows should match the length of
#' \code{new_data}, and the number of columns should match the dimensions of
#' the \code{xreg} argument in \code{object}. These rows extend
#' \code{object$xreg} and correspond to the periods immediately following it.
#' @param ... Other arguments are passed to \code{forecastfun}.
#'
#' @return
#' A refreshed object of class \code{"cpforecast"} with updated fields (e.g.,
#' \code{x}, \code{MEAN}, \code{ERROR}, \code{LOWER}, \code{UPPER}, and any
#' method-specific components), reflecting newly appended data and re-computed
#' cross-validation forecasts and conformal prediction intervals.
#'
#' @examples
#' # Simulate time series from an AR(2) model
#' library(forecast)
#' set.seed(1)
#' series <- arima.sim(n = 200, list(ar = c(0.8, -0.5)), sd = sqrt(1))
#'
#' # Cross-validation forecasting
#' far2 <- function(x, h, level) {
#'   Arima(x, order = c(2, 0, 0)) |>
#'     forecast(h = h, level)
#' }
#' fc <- cvforecast(series, forecastfun = far2, h = 3, level = 95,
#'                  forward = TRUE, initial = 1, window = 50)
#'
#' # Classical conformal prediction with equal weights
#' scpfc <- conformal(fc, method = "scp", symmetric = FALSE, ncal = 50, rolling = TRUE)
#'
#' # Update conformal prediction using newly available data
#' scpfc_update <- update(scpfc, forecastfun = far2, new_data = c(1.5, 0.8, 2.3))
#' print(scpfc_update)
#' summary(scpfc_update)
#'
#' @export
update.cpforecast <- function(
  object,
  new_data,
  forecastfun,
  new_xreg = NULL,
  ...
) {
  level <- object$level
  h <- dim(object$MEAN)[2]

  cvargs <- object$model$cvforecast$args
  if (is.null(cvargs)) {
    stop(
      "`object` carries no record of the `cvforecast` arguments; it was not ",
      "produced by scp(), acp(), pid() or acmcp()",
      call. = FALSE
    )
  }
  forward <- cvargs$forward
  window <- cvargs$window

  # Append new data
  if (!is.numeric(new_data) || !is.null(dim(new_data)) || !length(new_data)) {
    stop("`new_data` should be a non-empty numeric vector")
  }
  n_new <- length(new_data)
  x <- ts(
    c(object$x, new_data),
    start = start(object$x),
    frequency = frequency(object$x)
  )
  has_xreg <- "xreg" %in% names(object)
  if (has_xreg && is.null(new_xreg)) {
    stop("`new_xreg` is required because `object` contains `xreg`")
  }
  if (!has_xreg && !is.null(new_xreg)) {
    stop("`new_xreg` should only be supplied when `object` contains `xreg`")
  }
  if (has_xreg) {
    if (!is.numeric(new_xreg)) {
      stop("`new_xreg` should be a numeric matrix or vector")
    }
    new_xreg <- as.matrix(new_xreg)
    if (nrow(new_xreg) != n_new) {
      stop("`new_xreg` should have one row per observation in `new_data`")
    }
    if (ncol(new_xreg) != ncol(object$xreg)) {
      stop(
        "`new_xreg` and `object$xreg` should have the same number of columns"
      )
    }
    xreg <- ts(
      rbind(object$xreg, new_xreg),
      start = start(object$xreg),
      frequency = frequency(object$xreg)
    )
  } else {
    xreg <- NULL
  }

  # Model fitting and forecasting
  nfirst <- ifelse(forward, length(object$x) + 1L, length(object$x))
  nlast <- nfirst + n_new - 1L
  indx <- seq(nfirst, nlast, by = 1L)

  MEAN <- rbind(object$MEAN, matrix(NA, nrow = n_new, ncol = h)) |>
    ts(start = start(object$MEAN), frequency = frequency(object$MEAN))
  LOWER <- lapply(object$LOWER, function(lo) {
    rbind(lo, matrix(NA, nrow = n_new, ncol = h)) |>
      ts(start = start(lo), frequency = frequency(lo))
  })
  UPPER <- lapply(object$UPPER, function(up) {
    rbind(up, matrix(NA, nrow = n_new, ncol = h)) |>
      ts(start = start(up), frequency = frequency(up))
  })
  ERROR <- rbind(object$ERROR, matrix(NA, nrow = n_new, ncol = h)) |>
    ts(start = start(object$ERROR), frequency = frequency(object$ERROR))

  offset <- round((tsp(MEAN)[1L] - tsp(x)[1L]) * frequency(x))
  for (i in indx) {
    x_subset <- subset(
      x,
      start = ifelse(is.null(window), 1L, i - window + 1L),
      end = i
    )
    if (is.null(xreg)) {
      fc <- try(
        suppressWarnings(
          forecastfun(x_subset, h = h, level = level, ...)
        ),
        silent = TRUE
      )
    } else {
      xreg_subset <- subset(
        xreg,
        start = ifelse(is.null(window), 1L, i - window + 1L),
        end = i
      )
      xreg_future <- subset(
        xreg,
        start = i + 1L,
        end = i + h
      )
      fc <- try(
        suppressWarnings(
          forecastfun(
            x_subset,
            h = h,
            level = level,
            xreg = xreg_subset,
            newxreg = xreg_future,
            ...
          )
        ),
        silent = TRUE
      )
    }

    if (!is.element("try-error", class(fc))) {
      tm <- i - offset
      if (tm < 1L || tm + h > nrow(MEAN)) {
        stop(
          "cannot place the forecasts for origin ",
          i,
          " in `MEAN`",
          call. = FALSE
        )
      }
      MEAN[cbind(tm + 1:h, 1:h)] <- fc$mean
    }
  }
  ERROR[(nrow(ERROR) - n_new + 1):nrow(ERROR), ] <- new_data -
    MEAN[(nrow(ERROR) - n_new + 1):nrow(ERROR), ]

  # Update object info for conformal
  object$x <- x
  if (!is.null(xreg)) {
    object$xreg <- xreg
  }
  if (forward) {
    if (inherits(fc, "try-error")) {
      object$mean <- NULL
      warning(
        "the final (forward) model fit failed; forward forecasts are not available",
        call. = FALSE
      )
    } else {
      object$mean <- fc$mean
    }
  }
  object$MEAN <- MEAN
  object$ERROR <- ERROR
  object$LOWER <- LOWER
  object$UPPER <- UPPER
  object$forward <- forward
  if (object$method == "acp") {
    object$model$alpha_update <- lapply(
      object$model$alpha_update,
      function(alp) {
        lapply(alp, function(lv) {
          rbind(lv, matrix(NA, nrow = n_new, ncol = h)) |>
            ts(start = start(lv), frequency = frequency(lv))
        })
      }
    )
  }

  # Conformal prediction. Replay the resolved arguments, not the unevaluated
  # expressions in `object$call`.
  args <- object$model$args
  if (is.null(args)) {
    stop("`object` carries no record of its conformal arguments", call. = FALSE)
  }
  args$object <- object
  args$method <- object$method
  args$update <- TRUE
  out <- do.call(conformal, args)

  return(out)
}
