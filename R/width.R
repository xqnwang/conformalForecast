#' Calculate interval forecast width
#'
#' Calculate the mean width of prediction intervals on validation set (and test
#' set if applicable). If \code{window} is not \code{NULL}, the rolling mean
#' matrix for interval width is also returned. If \code{includemedian} is \code{TRUE},
#' the median interval width information will also be returned.
#'
#' @aliases print.width
#'
#' @param object An object of class \code{cvforecast} or \code{cpforecast}.
#' @param ... Additional inputs if \code{object} is missing.
#' @param level Target confidence level for prediction intervals.
#' @param includemedian If \code{TRUE}, the median interval width will also be returned.
#' @param window If not \code{NULL}, the rolling mean (and median if applicable)
#' matrix for interval width is also returned.
#' @param na.rm A logical indicating whether \code{NA} values should be stripped
#' before the mean or rolling mean computation proceeds.
#'
#' @return A list of class \code{"width"} with the following components:
#' \item{width}{Forecast interval width as a multivariate time series, where the \eqn{h}th
#' column holds the interval width for forecast horizon \eqn{h}. The time index
#' corresponds to the period for which the forecast is produced.}
#' \item{mean}{The mean interval width across the validation set (and test set
#' if applicable).}
#' \item{rollmean}{If \code{window} is not NULL, a rolling mean interval width
#' will be returned.}
#' \item{median}{The median interval width across the validation set (and test set
#' if applicable).}
#' \item{rollmedian}{If \code{window} is not NULL, a rolling median interval width
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
#' # Mean and rolling mean width for interval forecasts on validation and test sets
#' wid_fc <- width(fc, level = 95, window = 100)
#' str(wid_fc)
#'
#' @export
width <- function(object, ..., level = 95, includemedian = FALSE, window = NULL, na.rm = TRUE) {
  # Check inputs
  if (level > 0 && level < 1) {
    level <- 100 * level
  } else if (level < 0 || level > 99.99) {
    stop("confidence limit out of range")
  }
  dots <- rlang::dots_list(...)
  if (missing(object)) {
    if (any(!(c("LOWER", "UPPER") %in% names(dots))))
      stop("LOWER, and UPPER are required for interval width calculation")
  } else {
    if (any(!(c("LOWER", "UPPER") %in% names(object))))
      stop("LOWER, and UPPER are required for interval width calculation")
    if (!(level %in% object$level))
      stop("no interval forecasts of target confidence level in object")
    levelname <- paste0(level, "%")
    LOWER <- object$LOWER[[levelname]]
    UPPER <- object$UPPER[[levelname]]
  }
  lower <- LOWER
  upper <- UPPER
  horizon <- ncol(lower)
  period <- frequency(lower)

  # Match time
  tspl <- tsp(lower)
  tspu <- tsp(upper)
  start <- max(tspl[1], tspu[1])
  end <- min(tspl[2], tspu[2])

  lower <- window(lower, start = start, end = end)
  upper <- window(upper, start = start, end = end)
  n <- nrow(lower)

  # Width matrix
  widmat <- (upper- lower) |>
    ts(start = start, end = end, frequency = period)
  colnames(widmat) <- colnames(lower)

  out <- list(
    width = widmat
  )

  # Mean width
  out$mean <- apply(widmat, 2, mean, na.rm = na.rm)

  # Rolling mean width
  if (!is.null(window)) {
    if (window >= n)
      stop("the `window` argument should be smaller than the total period of interest")
    out$rollmean <- apply(widmat, 2, zoo::rollmean, k = window, na.rm = na.rm) |>
        ts(end = end, frequency = period)
  }

  # Median
  if (includemedian) {
    # Median width
    out$median <- apply(widmat, 2, median, na.rm = na.rm)

    # Rolling median width
    if (!is.null(window)) {
      if (window >= n)
        stop("the `window` argument should be smaller than the total period of interest")
      out$rollmedian <- apply(widmat, 2, zoo::rollmedian, k = window, na.rm = na.rm) |>
        ts(end = end, frequency = period)
    }
  }

  return(structure(out, class = "width"))
}

#' @export
print.width <- function(x, ...) {
  cat("Mean width:\n")
  print(x$mean)

  if ("median" %in% names(x)) {
    cat("\nMedian width:\n")
    print(x$median)
  }
}
