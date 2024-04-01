#' @rdname point_measures
#' @export
ME <- function(resid, na.rm = TRUE) {
  mean(resid, na.rm = na.rm)
}

#' @rdname point_measures
#' @export
MAE <- function(resid, na.rm = TRUE, ...){
  mean(abs(resid), na.rm = na.rm)
}

#' @rdname point_measures
#' @export
MSE <- function(resid, na.rm = TRUE, ...){
  mean(resid ^ 2, na.rm = na.rm)
}

#' @rdname point_measures
#' @export
RMSE <- function(resid, na.rm = TRUE, ...){
  sqrt(MSE(resid, na.rm = na.rm))
}

#' @rdname point_measures
#' @export
MPE <- function(resid, actual, na.rm = TRUE, ...){
  mean(resid / actual * 100, na.rm = na.rm)
}

#' @rdname point_measures
#' @export
MAPE <- function(resid, actual, na.rm = TRUE, ...){
  mean(abs(resid / actual * 100), na.rm = na.rm)
}

#' @rdname point_measures
#' @export
MASE <- function(resid, train, demean = FALSE, na.rm = TRUE,
                 period, d = period == 1, D = period > 1, ...){
  if (D > 0) { # seasonal differencing
    train <- diff(train, lag = period, differences = D)
  }
  if (d > 0) {
    train <- diff(train, differences = d)
  }
  if(demean){
    train <- train - mean(train, na.rm = na.rm)
  }
  scale <- mean(abs(train), na.rm = na.rm)
  mean(abs(resid / scale), na.rm = na.rm)
}

#' @rdname point_measures
#' @export
RMSSE <- function(resid, train, demean = FALSE, na.rm = TRUE,
                  period, d = period == 1, D = period > 1, ...){
  if (D > 0) { # seasonal differencing
    train <- diff(train, lag = period, differences = D)
  }
  if (d > 0) {
    train <- diff(train, differences = d)
  }
  if(demean){
    train <- train - mean(train, na.rm = na.rm)
  }
  scale <- mean(train^2, na.rm = na.rm)
  sqrt(mean(resid^2 / scale, na.rm = na.rm))
}

#' Point estimate accuracy measures
#'
#' @param resid A vector of residuals from either the validation or test data.
#' @param actual A vector of responses matching the forecasts.
#' @param train A vector of responses used to train the model.
#' @param period The seasonal period of the data.
#' @param d Should the response model include a first difference?
#' @param D Should the response model include a seasonal difference?
#' @param demean Should the response be demeaned (for MASE and RMSSE)?
#' @param na.rm If \code{TRUE}, remove the missing values before calculating the
#' accuracy measure.
#' @param ... Additional arguments for each measure.
#'
#' @export
point_measures <- list(ME = ME, MAE = MAE, MSE = MSE, RMSE = RMSE, MPE = MPE,
                       MAPE = MAPE, MASE = MASE, RMSSE = RMSSE)

#' @rdname interval_measures
#' @export
winkler_score <- function(lower, upper, actual, level = 95, na.rm = TRUE, ...){
  if (level > 0 && level < 1) {
    level <- 100 * level
  } else if (level < 0 || level > 99.99) {
    stop("confidence limit out of range")
  }
  alpha <- 1 - level/100
  score <- ifelse(
    actual < lower,
    (upper - lower) + (2 / alpha) * (lower - actual),
    ifelse(
      actual > upper,
      (upper - lower) + (2 / alpha) * (actual - upper),
      # else
      upper - lower)
  )
  mean(score, na.rm = na.rm)
}

#' @rdname interval_measures
#' @export
MSIS <- function(lower, upper, actual, train, level = 95,
                 period, d = period == 1, D = period > 1,
                 na.rm = TRUE, ...) {
  if (D > 0) { # seasonal differencing
    train <- diff(train, lag = period, differences = D)
  }
  if (d > 0) {
    train <- diff(train, differences = d)
  }
  scale <- mean(abs(train), na.rm = na.rm)
  score <- winkler_score(lower = lower, upper = upper, actual = actual,
                         level = level, na.rm = na.rm)
  mean(score / scale, na.rm = na.rm)
}


#' Interval estimate accuracy measures
#'
#' @param lower A vector of lower bound of interval forecasts from either the
#' validation or test data.
#' @param upper A vector of upper bound of interval forecasts from either the
#' validation or test data.
#' @param actual A vector of responses matching the forecasts.
#' @param train A vector of responses used to train the model.
#' @param level The level of the forecast interval.
#' @param period The seasonal period of the data.
#' @param d Should the response model include a first difference?
#' @param D Should the response model include a seasonal difference?
#' @param na.rm If \code{TRUE}, remove the missing values before calculating the
#' accuracy measure.
#' @param ... Additional arguments for each measure.
#'
#' @export
interval_measures <- list(Winkler = winkler_score, MSIS = MSIS)

#' Accuracy measures for a cross-validation model and a conformal prediction model
#'
#' Return range of summary measures of the out-of-sample forecast accuracy.
#' If \code{x} is given, the function also measures test set forecast accuracy.
#' If \code{x} is not given, the function only produces validation set forecast
#' accuracy.
#'
#' The measures calculated are:
#' \itemize{
#'   \item ME: Mean Error
#'   \item MAE: Mean Absolute Error
#'   \item MSE: Mean Squared Error
#'   \item RMSE: Root Mean Squared Error
#'   \item MPE: Mean Percentage Error
#'   \item MAPE: Mean Absolute Percentage Error
#'   \item MASE: Mean Absolute Scaled Error
#'   \item RMSSE: Root Mean Squared Scaled Error
#'   \item winkler_score: Winkler Score
#'   \item MSIS: Mean Scaled Interval Score
#' }
#'
#' @param object An object of class \code{cvforecast} or \code{cpforecast}.
#' @param x An optional numerical vector containing actual values of the same
#' length as object$mean.
#' @param CV If \code{TRUE}, the cross-validation forecast accuracy will be returned.
#' @param period The seasonal period of the data.
#' @param measures A list of accuracy measure functions to compute (such as \link{point_measures} or \link{interval_measures}).
#' @param byhorizon If \code{TRUE}, forecast accuracy will be calculated for each
#' individual forecast horizon \code{h} separately.
#' @param ... Additional arguments depending on the specific method.
#'
#' @return Matrix giving mean out-of-sample forecast accuracy measures.
#'
#' @seealso \code{\link{point_measures}}, \code{\link{interval_measures}}
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
#' # Out-of-sample forecast accuracy on validation set
#' accuracy(fc, measures = point_measures, byhorizon = TRUE)
#' accuracy(fc, measures = interval_measures, level = 95, byhorizon = TRUE)
#'
#' # Out-of-sample forecast accuracy on test set
#' accuracy(fc, x = c(1, 0.5, 0), measures = interval_measures,
#'          level = 95, byhorizon = TRUE)
#'
#' @importFrom stats frequency tsp window
#' @importFrom rlang dots_list
#' @export
accuracy.default <- function(object, x, CV = TRUE, period = NULL,
                             measures = interval_measures,
                             byhorizon = FALSE, ...) {
  if (!any(is.element(class(object), c("cvforecast", "cpforecast"))))
    stop(paste("no accuracy method found for an object of class",class(object)))
  if (!is.list(measures))
    stop("the `measures` argument must contain a list of accuracy measures")

  # Confidence level for interval measures
  dots <- rlang::dots_list(...)
  level <- ifelse(
    "level" %in% names(dots),
    ifelse(dots$level > 0 && dots$level < 1, 100*dots$level, dots$level),
    ifelse(95 %in% object$level, 95, object$level[1]))
  if (!(level %in% object$level))
    stop("the `level` argument should exist in object$level")
  # Change name for interval measures to include level
  names(measures) <- sapply(1:length(measures), function(i) {
    ifelse("level" %in% names(formals(measures[[i]])),
           paste0(names(measures)[i], "_", level),
           names(measures)[i])
  })

  cvset <- (is.list(object))
  testset <- (!missing(x))
  if (!missing(x)) {
    xx <- x
  }
  if (!cvset && !testset)
    stop("unable to compute forecast accuracy measures")
  if (!CV)
    cvset <- FALSE

  if (is.null(period)) {
    if (testset) {
      period <- frequency(xx)
    }
    if (cvset) {
      period <- frequency(object$x)
    }
  } else {
    if (period != frequency(object$x)) {
      stop("the `period` argument does not match the frequency of object$x")
    }
  }

  if (cvset) {
    x <- object$x
    e <- object$ERROR

    lb <- object$LOWER[[paste0(level, "%")]]
    ub <- object$UPPER[[paste0(level, "%")]]

    horizon <- ncol(e)

    if (all(is.ts(x), is.ts(e), is.ts(lb), is.ts(ub))) {
      tspo <- tsp(x)
      tspe <- tsp(e)
      if (tspo[1] >= tspe[1])
        stop("the start time of object$ERROR should be later than that of object$x")
      end <- min(tspe[2], tspo[2])
      x <- window(x, start = tspo[1], end = end)
      e <- window(e, start = tspe[1], end = end) |> lagmatrix(-(0:(horizon-1)))
      lb <- window(lb, start = tspe[1], end = end) |> lagmatrix(-(0:(horizon-1)))
      ub <- window(ub, start = tspe[1], end = end) |> lagmatrix(-(0:(horizon-1)))
    } else {
      stop("both x and ERROR in object should be time-series objects")
    }
    n <- length(x) - nrow(e)

    if (byhorizon) {
      input <- 1:horizon
      cvrnames <- paste0("CV h=", 1:horizon)
    } else {
      input <- list(1:horizon)
      cvrnames <- paste0("CV")
    }

    cvout <- lapply(input, function(h) {
      out_h <- lapply(1:nrow(e), function(t) {
        all_args <- list(
          resid = e[t, h], actual = x[n+t-1+h], train = x[1:(n+t-1)], period = period,
          lower = lb[t, h], upper = ub[t, h], level = level, ...
        )
        sapply(measures, function(f) {
          f_args <- all_args[intersect(names(formals(f)), names(all_args))]
          tryCatch(
            {do.call(f, f_args)},
            error = function(e) return(NA_real_)
          )
        })
      })
      out_h <- do.call(rbind, out_h)
      apply(out_h, 2, mean, na.rm = TRUE)
    })
    cvout <- do.call(rbind, cvout)
    rownames(cvout) <- cvrnames
  } else {
    cvout <- NULL
  }

  if (testset) {
    x <- object$x
    ff <- object$mean
    if (length(ff) != length(xx))
      stop("the length of object$mean and x do not match")
    horizon <- length(xx)
    ee <- xx - ff

    lt <- object$lower[, paste0(level, "%")]
    ut <- object$upper[, paste0(level, "%")]

    if (byhorizon) {
      input <- 1:horizon
      testrnames <- paste0("Test h=", 1:horizon)
    } else {
      input <- list(1:horizon)
      testrnames <- paste0("Test")
    }

    testout <- lapply(input, function(h) {
      all_args <- list(
        resid = ee[h], actual = xx[h], train = x, period = period,
        lower = lt[h], upper = ut[h], level = level, ...
      )
      sapply(measures, function(f) {
        f_args <- all_args[intersect(names(formals(f)), names(all_args))]
        tryCatch(
          {do.call(f, f_args)},
          error = function(e) return(NA_real_)
        )
      })
    })
    testout <- do.call(rbind, testout)
    rownames(testout) <- testrnames
  } else {
    testout <- NULL
  }

  out <- rbind(cvout, testout)
  return(out)
}
