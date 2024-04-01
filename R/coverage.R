#' Calculate interval forecast coverage
#'
#' Calculate the mean coverage and the ifinn matrix for prediction intervals on
#' validation set (and test set if applicable). If \code{window} is not \code{NULL},
#' the rolling mean matrix for coverage is also returned.
#'
#' @aliases print.coverage
#'
#' @param object An object of class \code{cvforecast} or \code{cpforecast}.
#' @param ... Additional inputs if \code{object} is missing.
#' @param level Target confidence level for prediction intervals.
#' @param window If not \code{NULL}, the rolling mean matrix for coverage is also returned.
#' @param na.rm A logical indicating whether \code{NA} values should be stripped
#' before the mean or rolling mean computation proceeds.
#'
#' @return A list of class \code{"coverage"} with the following components:
#' \item{mean}{The mean coverage across the validation set (and test set if applicable).}
#' \item{ifinn}{A indicator matrix as a multivariate time series, where the \eqn{h}th column
#' holds the coverage for forecast horizon \eqn{h}. The time index
#' corresponds to the period for which the forecast is produced.}
#' \item{rollmean}{If \code{window} is not NULL, a rolling mean matrix for coverage
#' will be returned.}
#'
#' @examples
#' # Simulate time series from an AR(2) model
#' library(forecast)
#' series <- arima.sim(n = 1000, list(ar = c(0.8, -0.5)), sd = sqrt(1))
#' series <- as.numeric(series)
#'
#' # Cross-validation forecasting with a rolling window of length 100
#' far2 <- function(x, h, level) {
#'   Arima(x, order = c(2, 0, 0)) |>
#'     forecast(h = h, level)
#' }
#' fc <- cvforecast(series, forecastfun = far2, h = 3, level = c(80, 95),
#'                  forward = TRUE, initial = 1, window = 100)
#'
#' # Mean and rolling mean coverage for interval forecasts on validation and test sets
#' cov_fc <- coverage(fc, level = 95, window = 100)
#' str(cov_fc)
#'
#' @importFrom stats window
#' @importFrom zoo rollmean
#' @export
coverage <- function(object, ..., level = 95, window = NULL, na.rm = TRUE) {
  # Check inputs
  if (level > 0 && level < 1) {
    level <- 100 * level
  } else if (level < 0 || level > 99.99) {
    stop("confidence limit out of range")
  }
  dots <- rlang::dots_list(...)
  if (missing(object)) {
    if (any(!(c("x", "LOWER", "UPPER") %in% names(dots))))
      stop("x, LOWER, and UPPER are required for coverage calculation")
  } else {
    if (any(!(c("x", "LOWER", "UPPER") %in% names(object))))
      stop("x, LOWER, and UPPER are required for coverage calculation")
    if (!(level %in% object$level))
      stop("no interval forecasts of target confidence level in object")
    levelname <- paste0(level, "%")
    x <- object$x
    LOWER <- object$LOWER[[levelname]]
    UPPER <- object$UPPER[[levelname]]
  }
  lower <- LOWER
  upper <- UPPER
  horizon <- ncol(lower)
  period <- frequency(object$x)
  x <- ts(matrix(rep(object$x, horizon), ncol = horizon, byrow = FALSE),
          start = start(object$x),
          frequency = period)

  # Match time
  tspx <- tsp(x)
  tspl <- tsp(lower)
  tspu <- tsp(upper)
  start <- max(tspx[1], tspl[1], tspu[1])
  end <- min(tspx[2], tspl[2], tspu[2])

  x <- window(x, start = start, end = end)
  lower <- window(lower, start = start, end = end)
  upper <- window(upper, start = start, end = end)
  n <- nrow(x)

  # If coverage matrix
  covmat <- (lower <= x & x <= upper) |>
    ts(start = start, end = end, frequency = period)
  colnames(covmat) <- colnames(lower)

  # Mean coverage
  covmean <- apply(covmat, 2, mean, na.rm = na.rm)

  # Rolling mean coverage
  if (!is.null(window)) {
    if (window >= n)
      stop("the `window` argument should be smaller than the total period of interest")
    covrmean <- apply(covmat, 2, zoo::rollmean, k = window, na.rm = na.rm) |>
      ts(end = end, frequency = period)
  }

  out <- list(
    mean = covmean,
    ifinn = covmat
  )
  if (!is.null(window)) out <- append(out, list(rollmean = covrmean))
  return(structure(out, class = "coverage"))
}

#' @export
print.coverage <- function(x, ...) {
  print(x$mean)
}
