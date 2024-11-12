#' Conformal prediction
#'
#' This function allows you to specify the method used to perform conformal prediction.
#'
#' @param object An object of class \code{"cvforecast"}. It must have an argument
#' \code{x} for original univariate time series, an argument \code{MEAN} for
#' point forecasts and \code{ERROR} for forecast errors on validation set.
#' See the results of a call to \code{\link{cvforecast}}.
#' @param method A character string specifying the conformal method to be applied.
#' Possible options include \code{"scp"} (\link{scp}), \code{"acp"} (\link{acp}), \code{"pid"} (\link{pid}), and \code{"mcp"} (\link{mcp}).
#' @param ... Additional arguments to be passed to the selected conformal method.
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
#' summary(scpfc)
#'
#' # ACP with asymmetric nonconformity scores and rolling calibration sets
#' acpfc <- conformal(fc, method = "acp", symmetric = FALSE, gamma = 0.005, ncal = 100, rolling = TRUE)
#' summary(acpfc)
#'
#' @export
conformal <- function(object, ...) {
  UseMethod("conformal")
}

#' @rdname conformal
#' @export
conformal.cvforecast <- function(object, method = c("scp", "acp", "pid", "mcp"), ...) {
  # Match the method argument
  method <- match.arg(method)

  # Call the respective function based on the method
  result <- switch(method,
                   scp = scp(object, ...),
                   acp = acp(object, ...),
                   pid = pid(object, ...),
                   mcp = mcp(object, ...))

  # Return the result of the chosen method
  return(result)
}
