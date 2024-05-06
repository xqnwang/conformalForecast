#' Time series cross-validation forecasting
#'
#' Compute forecasts and other information by applying
#' \code{forecastfun} to subsets of the time series \code{y} using a
#' rolling forecast origin.
#'
#' Let \code{y} denote the time series \eqn{y_1,\dots,y_T}{y[1:T]} and
#' let \eqn{t_0} denote the initial period.
#'
#' Suppose \code{forward = TRUE}. If \code{window} is \code{NULL},
#' \code{forecastfun} is applied successively to the subset time series
#' \eqn{y_{1},\dots,y_t}{y[1:t]}, for \eqn{t=t_0,\dots,T},
#' generating forecasts \eqn{\hat{y}_{t+1|t},\dots,\hat{y}_{t+h|t}}{f[t+1:h]}. If \code{window} is not
#' \code{NULL} and has a length of \eqn{t_w}, then \code{forecastfun} is applied
#' successively to the subset time series \eqn{y_{t-t_w+1},\dots,y_{t}}{y[(t-t_w+1):t)]},
#' for \eqn{t=\max(t_0, t_w),\dots,T}.
#'
#' If \code{forward} is \code{FALSE}, the last observation used for training will
#' be \eqn{y_{T-1}}.
#'
#' @aliases print.cvforecast summary.cvforecast print.summary.cvforecast
#'
#' @param y Univariate time series.
#' @param forecastfun Function to return an object of class \code{"forecast"}.
#' Its first argument must be a univariate time series, and it must have an
#' argument \code{h} for the forecast horizon and an argument \code{level} for
#' the confidence level for prediction intervals. If exogenous predictors are used,
#' then it must also have \code{xreg} and \code{newxreg} arguments corresponding
#' to the training and test periods, respectively.
#' @param h Forecast horizon.
#' @param level Confidence level for prediction intervals.
#' @param forward If \code{TRUE}, the final forecast origin for forecasting is
#' \eqn{y_T}. Otherwise, the final forecast origin is \eqn{y_{T-1}}.
#' @param xreg Exogenous predictor variables passed to \code{forecastfun} if required.
#' It should be of the same size as \code{y}+\code{forward}*\code{h}, otherwise,
#' \code{NA} padding or subsetting will be applied.
#' @param initial Initial period of the time series where no cross-validation
#' forecasting is performed.
#' @param window Length of the rolling window. If \code{NULL}, a rolling window
#' will not be used.
#' @param ... Other arguments are passed to \code{forecastfun}.
#'
#' @return A list of class \code{c("cvforecast", "forecast")} with components:
#' \item{x}{The original time series.}
#' \item{series}{The name of the series \code{x}.}
#' \item{method}{A character string "cvforecast".}
#' \item{fit_times}{The number of times the model is fitted in cross-validation.}
#' \item{MEAN}{Point forecasts as a multivariate time series, where the \eqn{h}th column
#' holds the point forecasts for forecast horizon \eqn{h}. The time index
#' corresponds to the period for which the forecast is produced.}
#' \item{ERROR}{Forecast errors given by
#' \eqn{e_{t+h|t} = y_{t+h}-\hat{y}_{t+h|t}}{e[t+h] = y[t+h]-f[t+h]}.}
#' \item{LOWER}{A list containing lower bounds for prediction intervals for
#' each \code{level}. Each element within the list will be a multivariate time
#' series with the same dimensional characteristics as \code{MEAN}.}
#' \item{UPPER}{A list containing upper bounds for prediction intervals for
#' each \code{level}. Each element within the list will be a multivariate time
#' series with the same dimensional characteristics as \code{MEAN}.}
#' \item{level}{The confidence values associated with the prediction intervals.}
#' \item{call}{The matched call.}
#' If \code{forward} is \code{TRUE}, the components \code{mean}, \code{lower},
#' \code{upper}, and \code{model} will also be returned, showing the information
#' about the final fitted model and forecasts using all available observations, see
#' e.g. \code{\link[forecast]{forecast.ets}} for more details.
#'
#' @examples
#' # Simulate time series from an AR(2) model
#' library(forecast)
#' series <- arima.sim(n = 1000, list(ar = c(0.8, -0.5)), sd = sqrt(1))
#'
#' # Example with a rolling window of length 100
#' far2 <- function(x, h, level) {
#'   Arima(x, order = c(2, 0, 0)) |>
#'     forecast(h = h, level)
#' }
#' fc <- cvforecast(series, forecastfun = far2, h = 3, level = c(80, 95),
#'                  forward = TRUE, initial = 1, window = 100)
#' print(fc)
#' summary(fc)
#'
#' # Example with exogenous predictors
#' far2_xreg <- function(x, h, level, xreg, newxreg) {
#'   Arima(x, order=c(2, 0, 0), xreg = xreg) |>
#'     forecast(h = h, level = level, xreg = newxreg)
#' }
#' fc_xreg <- cvforecast(series, forecastfun = far2_xreg, h = 3, level = c(80, 95),
#'                       forward = TRUE, xreg = matrix(rnorm(2006), ncol = 2, nrow = 1003),
#'                       initial = 1, window = 100)
#'
#' @export
cvforecast <- function(y, forecastfun, h = 1, level = c(80, 95),
                       forward = TRUE, xreg = NULL, initial = 1, window = NULL, ...) {
  # Check whether there are non-existent arguments
  all.args <- names(formals())
  user.args <- names(match.call())[-1L] # including arguments passed to 3 dots
  check <- user.args %in% all.args
  if (!all(check)) {
    error.args <- user.args[!check]
    warning(sprintf("The non-existent %s arguments will be ignored.", error.args))
  }

  # Check input univariate time series
  if (any(class(y) %in% c("data.frame", "list", "matrix", "mts"))) {
    stop("y should be a univariate time series")
  }
  seriesname <- deparse(substitute(y))
  y <- as.ts(y)
  n <- length(y)

  # Check confidence level
  if (min(level) > 0 && max(level) < 1) {
    level <- 100 * level
  } else if (min(level) < 0 || max(level) > 99.99) {
    stop("confidence limit out of range")
  }
  level <- sort(level)

  # Check other inputs
  if (h <= 0)
    stop("forecast horizon out of bounds")
  if (initial < 1 | initial > n)
    stop("initial period out of bounds")
  if (initial == n && !forward)
    stop("initial period out of bounds")
  if (!is.null(window)) {
    if (window < 1 | window > n)
      stop("window out of bounds")
    if (window == n && !forward)
      stop("window out of bounds")
  }
  if (!is.null(xreg)) {
    if(!is.numeric(xreg))
      stop("xreg should be a numeric matrix or a numeric vector")
    xreg <- ts(as.matrix(xreg),
               start = start(y),
               frequency = frequency(y))

    if (nrow(xreg) < n)
      stop("xreg should be at least of the same size as y")
    if (nrow(xreg) < n + forward * h)
      # Pad xreg with NAs
      xreg <- ts(rbind(xreg, matrix(NA, nrow=n+forward*h-nrow(xreg), ncol=ncol(xreg))),
                 start = start(y),
                 frequency = frequency(y))
    if (nrow(xreg) > n + forward * h) {
      warning(sprintf("only first %s rows of xreg are being used", n + forward * h))
      xreg <- subset(xreg, start = 1L, end = n + forward * h)
    }
    if (is.null(colnames(xreg))) {
      colnames(xreg) <- if (ncol(xreg) == 1) "xreg" else paste0("xreg", 1:ncol(xreg))
    }
  }

  N <- ifelse(forward, n + h, n + h - 1L)
  nlast <- ifelse(forward, n, n - 1L)
  nfirst <- ifelse(is.null(window), initial, max(window, initial))
  indx <- seq(nfirst, nlast, by = 1L)
  fit_times <- length(indx)

  pf <- err <- `colnames<-` (
    ts(matrix(NA_real_, nrow = N, ncol = h),
       start = start(y),
       frequency = frequency(y)),
    paste0("h=", 1:h))
  lower <- upper <- `names<-` (rep(list(pf), length(level)), paste0(level, "%"))
  out <- list(
    x = y
  )

  for (i in indx) {
    y_subset <- subset(
      y,
      start = ifelse(is.null(window), 1L, i - window + 1L),
      end = i)

    if (is.null(xreg)) {
      fc <- try(suppressWarnings(
        forecastfun(y_subset, h = h, level = level, ...)
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
        forecastfun(y_subset, h = h, level = level,
                    xreg = xreg_subset, newxreg = xreg_future, ...)
        ), silent = TRUE)
    }

    if (!is.element("try-error", class(fc))) {
      pf[i,] <- fc$mean
      err[i,] <- y[i + 1:h] - fc$mean
      for (l in level) {
        levelname <- paste0(l, "%")
        lower[[levelname]][i,] <- fc$lower[, levelname]
        upper[[levelname]][i,] <- fc$upper[, levelname]
      }
    }
  }

  out$series <- seriesname
  out$method <- paste("cvforecast")
  out$fit_times <- fit_times
  out$MEAN <- lagmatrix(pf, 1:h) |> window(start = time(pf)[nfirst + 1L])
  out$ERROR <- lagmatrix(err, 1:h) |> window(start = time(err)[nfirst + 1L], end = time(err)[n])
  out$LOWER <- lapply(lower,
                      function(low) lagmatrix(low, 1:h) |>
                        window(start = time(low)[nfirst + 1L]))
  out$UPPER <- lapply(upper,
                      function(up) lagmatrix(up, 1:h) |>
                        window(start = time(up)[nfirst + 1L]))
  out$level <- level
  out$call <- match.call()
  # The final forecasting model output if forward is TRUE
  if (forward) {
    out$mean <- fc$mean
    out$lower <- fc$lower
    out$upper <- fc$upper
    out$model <- fc$model
  }

  return(structure(out, class = c("cvforecast", "forecast")))
}

#' @export
print.cvforecast <- function(x, ...) {
  cat(paste("Cross-validation\n\n"))
  if (!is.null(x$call)) {
    cat(paste("Call:\n"))
    for (i in 1:length(deparse(x$call))) {
      cat(paste("", deparse(x$call)[i]), "\n")
    }
    cat(paste("\n"))
  }

  cat(paste("", "fit_times =", x$fit_times,
            ifelse("model" %in% names(x), "(the forward step included)", ""), "\n"))

  if ("model" %in% names(x)) {
    cat(paste("\nForecasts of the forward step:\n"))
    NextMethod()
  }
}

#' @export
summary.cvforecast <- function(object, ...) {
  class(object) <- c("summary.cvforecast", class(object))
  object
}

#' @export
print.summary.cvforecast <- function(x, ...) {
  NextMethod()
  cat("\nCross-validation error measures:\n")
  print(round(
    accuracy.default(x, measures = c(point_measures, interval_measures),
                     byhorizon = FALSE),
    digits = 3))
}
