#' Update and repeform cross-validation forecasting and conformal prediction
#'
#' Update conformal prediction intervals and other information by applying the
#' \code{cvforecast} and \code{conformal} functions.
#'
#' @param object An object of class \code{"cpforecast"}.
#' @param new_data A vector of newly available data.
#' @param forecastfun Function to return an object of class \code{"forecast"}.
#' Its first argument must be a univariate time series, and it must have an
#' argument \code{h} for the forecast horizon and an argument \code{level} for
#' the confidence level for prediction intervals. If exogenous predictors are used,
#' then it must also have \code{xreg} and \code{newxreg} arguments corresponding
#' to the training and test periods, respectively.
#' @param new_xreg Newly available exogenous predictor variables passed to
#' \code{forecastfun} if required. The number of rows should match the length of
#' \code{new_data}, and the number of columns should match the dimensions of
#' the \code{xreg} argument in \code{object}.
#' @param ... Other arguments are passed to \code{forecastfun}.
#'
#' @examples
#' # Simulate time series from an AR(2) model
#' library(forecast)
#' series <- arima.sim(n = 1000, list(ar = c(0.8, -0.5)), sd = sqrt(1))
#'
#' # Cross-validation forecasting
#' far2 <- function(x, h, level) {
#'   Arima(x, order = c(2, 0, 0)) |>
#'     forecast(h = h, level)
#' }
#' fc <- cvforecast(series, forecastfun = far2, h = 3, level = c(80, 95),
#'                  forward = TRUE, initial = 1, window = 100)
#'
#' # Classical conformal prediction with equal weights
#' scpfc <- conformal(fc, method = "scp", symmetric = FALSE, ncal = 100, rolling = TRUE)
#'
#' # Update conformal prediction using newly available data
#' scpfc_update <- update(scpfc, forecastfun = far2, new_data = c(1.5, 0.8, 2.3))
#' print(scpfc_update)
#' summary(scpfc_update)
#'
#' @export
update.cpforecast <- function(object, new_data, forecastfun, new_xreg = NULL, ...) {
  level <- object$level
  alpha <- 1 - level / 100
  h <- dim(object$MEAN)[2]
  forward <- "mean" %in% names(object)

  # Append new data
  n_new <- length(new_data)
  x <- ts(c(object$x, new_data), start = start(object$x), frequency = frequency(object$x))
  if (!is.null(new_xreg) & ("xreg" %in% names(object))) {
    if(!is.numeric(new_xreg))
      stop("'new_xreg' should be a numeric matrix or a numeric vector")
    new_xreg <- as.matrix(new_xreg)
    if (nrow(new_xreg) != n_new)
      stop("the size of 'new_xreg' must match that of 'new_data'")
    xreg <- ts(rbind(object$xreg, new_xreg),
               start = start(object$xreg),
               frequency = frequency(object$xreg))
  } else {
    xreg <- NULL
  }

  # Info required for model fitting
  cvcall <- object$model$cvforecast$call
  initial <- eval(get_call_arg_with_defaults(cvforecast, cvcall, "initial"), .GlobalEnv)
  window <- eval(get_call_arg_with_defaults(cvforecast, cvcall, "window"), .GlobalEnv)

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

  for (i in indx) {
    x_subset <- subset(
      x,
      start = ifelse(is.null(window), 1L, i - window + 1L),
      end = i)
    if (is.null(xreg)) {
      fc <- try(suppressWarnings(
        forecastfun(x_subset, h = h, level = level, ...)
      ), silent = TRUE)
    } else {
      xreg_subset <- subset(
        xreg,
        start = ifelse(is.null(window), 1L, i - window + 1L),
        end = i)
      xreg_future <- subset(
        xreg,
        start = i + 1L,
        end = i + h)
      fc <- try(suppressWarnings(
        forecastfun(x_subset, h = h, level = level,
                    xreg = xreg_subset, newxreg = xreg_future, ...)
      ), silent = TRUE)
    }

    if (!is.element("try-error", class(fc))) {
      tm <- which(tail(as.numeric(time(x_subset)), 1) == as.numeric(time(MEAN)))
      MEAN[cbind(tm + 1:h, 1:h)] <- fc$mean
    }
  }
  ERROR[(nrow(ERROR)-n_new+1):nrow(ERROR),] <- new_data - MEAN[(nrow(ERROR)-n_new+1):nrow(ERROR),]

  # Update object info for conformal
  object$x <- x
  if (!is.null(xreg)) object$xreg <- xreg
  if (forward) object$mean <- fc$mean
  object$MEAN <- MEAN
  object$ERROR <- ERROR
  object$LOWER <- LOWER
  object$UPPER <- UPPER
  object$forward <- forward
  if (object$method == "acp") {
    object$model$alpha_update <- lapply(object$model$alpha_update, function(alp){
      lapply(alp, function(lv){
        rbind(lv, matrix(NA, nrow = n_new, ncol = h)) |>
          ts(start = start(lv), frequency = frequency(lv))
      })
    })
  }

  # Conformal prediction
  args <- as.list(object$call)[-1]
  args$object <- object
  args$method <- object$method
  args$update <- TRUE
  out <- do.call(conformal, args)

  return(out)
}
