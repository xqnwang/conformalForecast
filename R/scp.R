#' Classical split conformal prediction method
#'
#' Compute prediction intervals and other information obtained by
#' applying the classical split conformal prediction method.
#'
#' Consider a vector \eqn{s_{t+h|t}} that contains the conformity scores for the
#' \eqn{h}-step-ahead forecasts.
#'
#' If \code{symmetric} is \code{TRUE}, \eqn{s_{t+h|t}=|e_{t+h|t}|}.
#' When \code{rolling} is \code{FALSE}, the \eqn{\alpha}-quantile
#' \eqn{\hat{q}_{t+h|t}} are computed successively on expanding calibration sets
#' \eqn{s_{1+h|1},\dots,s_{t|t-h}}, for \eqn{t=\mathrm{ncal}+h-1,\dots,T}. Then the
#' prediction intervals will be
#' \eqn{[\hat{y}_{t+h|t}-\hat{q}_{t+h|t}, \hat{y}_{t+h|t}+\hat{q}_{t+h|t}]}.
#' When \code{rolling} is \code{TRUE}, the calibration sets will be of same length
#' \code{ncal}.
#'
#' If \code{symmetric} is \code{FALSE}, \eqn{s_{t+h|t}^{u}=e_{t+h|t}} for upper
#' interval bounds and \eqn{s_{t+h|t}^{l} = -e_{t+h|t}} for lower bounds.
#' Instead of computing \eqn{\alpha}-quantile, \eqn{\alpha/2}-quantiles for lower
#' bound (\eqn{\hat{q}_{t+h|t}^{l}}) and upper bound (\eqn{\hat{q}_{t+h|t}^{u}})
#' are calculated based on their conformity scores, respectively.
#' Then the prediction intervals will be
#' \eqn{[\hat{y}_{t+h|t}-\hat{q}_{t+h|t}^{l}, \hat{y}_{t+h|t}+\hat{q}_{t+h|t}^{u}]}.
#'
#' @aliases print.scp summary.scp print.summary.scp
#'
#' @param object An object of class "\code{cvforecast}". It must have an argument
#' \code{x} for original univariate time series, an argument \code{MEAN} for
#' point forecasts and \code{ERROR} for forecast errors. See the results of a call
#' to \code{\link{cvforecast}}.
#' @param alpha A numeric vector of significance levels to achieve a desired
#' coverage level \eqn{1-\alpha}.
#' @param symmetric If \code{TRUE}, symmetric conformity scores (i.e. \eqn{|e_{t+h|t}|})
#' are used. If \code{FALSE}, asymmetric conformity scores (i.e. \eqn{e_{t+h|t}})
#' are used, and then upper bounds and lower bounds are produced separately.
#' @param ncal Length of the calibration set. If \code{rolling = FALSE}, it denotes
#' the initial period of calibration sets. Otherwise, it indicates
#' the period of every rolling calibration set.
#' @param rolling If \code{TRUE}, a rolling window strategy will be adopted to
#' form the calibration set. Otherwise, expanding window strategy will be used.
#' @param quantiletype An integer between 1 and 9 determining the type of
#' quantile estimator to be used. Types 1 to 3 are for discontinuous quantiles,
#' types 4 to 9 are for continuous quantiles. See the
#' \code{\link[ggdist]{weighted_quantile}} function in the ggdist package.
#' @param weightfun Function to return a vector of weights used for sample quantile
#' computation. Its first argument must be an integer indicating the number of
#' observations for which weights are generated. If \code{NULL}, equal weights
#' will be used for sample quantile computation. Currently, only non-data-dependent
#' weights are supported.
#' @param kess If \code{TRUE}, Kish's effective sample size is used for sample
#' quantile computation.
#' @param na.rm If \code{TRUE}, corresponding entries in sample values and weights
#' are removed if either is \code{NA} when calculating sample quantile.
#' @param ... Other arguments are passed to \code{weightfun}.
#'
#' @return A list of class \code{c("scp", "cvforecast", "forecast")}
#' with the following components:
#' \item{x}{The original time series.}
#' \item{series}{The name of the series \code{x}.}
#' \item{method}{A character string "scp".}
#' \item{cp_times}{The number of times the conformal prediction is performed in
#' cross-validation.}
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
#' \item{model}{A list containing information abouth the conformal prediction model.}
#' If \code{mean} is included in the \code{object}, the components \code{mean}, \code{lower},
#' \code{upper}, and \code{model} will also be returned, showing the information
#' about the forecasts generated using all available observations.
#'
#' @seealso \code{\link[ggdist]{weighted_quantile}}
#'
#' @examples
#' # Simulate time series from an AR(2) model
#' library(forecast)
#' series <- arima.sim(n = 1000, list(ar = c(0.8, -0.5)), sd = sqrt(1))
#' series <- as.numeric(series)
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
#' scpfc <- scp(fc, symmetric = FALSE, ncal = 100, rolling = TRUE)
#' print(scpfc)
#' summary(scpfc)
#'
#' # Classical conformal prediction with exponential weights
#' expweight <- function(n) {
#'   0.99^{n+1-(1:n)}
#' }
#' scpfc_exp <- scp(fc, symmetric = FALSE, ncal = 100, rolling = TRUE,
#'                  weightfun = expweight, kess = TRUE)
#'
#' @importFrom ggdist weighted_quantile
#' @export
scp <- function(object, alpha = 1 - 0.01 * object$level,
                symmetric = FALSE, ncal = 10, rolling = FALSE,
                quantiletype = 1, weightfun = NULL, kess = FALSE, na.rm = TRUE,
                ...) {
  # Check inputs
  if (any(alpha >= 1 | alpha <= 0))
    stop("`alpha` should be in (0, 1)")
  if (ncal < 10)
    stop("length of calibration period, `ncal`, should at least be 10")
  if (!quantiletype %in% 1:9)
    stop("`quantiletype` is invalid. It must be in 1:9.")
  if (is.null(weightfun)) {
    # Equal weights
    weightfun <- function(n) rep(1, n)
  }
  if (kess) {
    # Kish's effective sample size for sample quantile computation
    kess <- function(w) sum(w)^2 / sum(w^2)
  } else {
    kess <- NULL
  }

  alpha <- sort(alpha, decreasing = TRUE)
  level <- 100 * (1 - alpha)
  pf <- ts(as.matrix(object$MEAN),
           start = start(object$MEAN),
           frequency = frequency(object$MEAN))
  errors <- ts(as.matrix(object$ERROR),
               start = start(object$ERROR),
               frequency = frequency(object$ERROR))
  horizon <- ncol(pf)
  n <- nrow(pf)

  if (ncal > nrow(errors))
    stop("`ncal` is larger than the number of rows in object$ERROR")

  namatrix <- ts(matrix(NA_real_, nrow = n, ncol = horizon),
                 start = start(pf),
                 frequency = frequency(pf))
  colnames(namatrix) <- paste0("h=", seq(horizon))
  lower <- upper <- `names<-` (rep(list(namatrix), length(alpha)),
                               paste0(level, "%"))

  out <- list(
    x = object$x,
    series = object$series
  )

  for (h in seq(horizon)) {
    indx <- seq(ncal+h-1, nrow(errors), by = 1L)

    for (t in indx) {
      errors_subset <- subset(
        errors[, h],
        start = ifelse(!rolling, 1, t - ncal + 1L),
        end = t)

      weight_subset <- weightfun(length(errors_subset) + 1L, ...)

      if (symmetric) {
        q_lo <- q_up <- ggdist::weighted_quantile(
          x = abs(c(errors_subset, Inf)),
          probs = 1 - alpha,
          weights = weight_subset,
          n = kess,
          type = quantiletype,
          na.rm = na.rm)
      } else {
        q_lo <- ggdist::weighted_quantile(
          x = - c(errors_subset, Inf),
          probs = 1 - alpha/2,
          weights = weight_subset,
          n = kess,
          type = quantiletype,
          na.rm = na.rm)
        q_up <- ggdist::weighted_quantile(
          x = c(errors_subset, Inf),
          probs = 1 - alpha/2,
          weights = weight_subset,
          n = kess,
          type = quantiletype,
          na.rm = na.rm)
      }
      for (i in seq(length(alpha))) {
        lower[[i]][t+h, h] <- pf[t+h, h] - q_lo[i]
        upper[[i]][t+h, h] <- pf[t+h, h] + q_up[i]
      }
    }
  }

  out$method <- paste("scp")
  out$cp_times <- length(indx)
  out$MEAN <- object$MEAN
  out$ERROR <- object$ERROR
  out$LOWER <- lower
  out$UPPER <- upper
  out$level <- level
  out$call <- match.call()
  if ("mean" %in% names(object)) {
    out$mean <- object$mean
    out$lower <- extract_final(lower, nrow = n, ncol = horizon, bench = out$mean)
    out$upper <- extract_final(upper, nrow = n, ncol = horizon, bench = out$mean)
  }
  out$model$method <- out$method
  out$model$call <- match.call()
  out$model$alpha <- alpha
  out$model$symmetric <- symmetric

  return(structure(out, class = c("scp", "cpforecast", "forecast")))
}

# Extract final step forecasts from x and copy attributes from bench to it
extract_final <- function(x, nrow, ncol, bench) {
  x <- sapply(
    names(x),
    function(xx)
      sapply(ncol:1 - 1, function(h) as.numeric(x[[xx]][nrow-h, ncol-h]))
    , simplify = FALSE)
  x <- do.call(cbind, x)
  copy_msts(bench, x)
}

## Copied from forecast:::copy_msts
copy_msts <- function(x, y) {
  if(NROW(x) > NROW(y)) {
    # Pad y with initial NAs
    if(NCOL(y) == 1) {
      y <- c(rep(NA, NROW(x) - NROW(y)), y)
    } else {
      y <- rbind(matrix(NA, ncol=NCOL(y), nrow = NROW(x) - NROW(y)), y)
    }
  } else if(NROW(x) != NROW(y)) {
    stop("x and y should have the same number of observations")
  }
  if(NCOL(y) > 1) {
    class(y) <- c("mts", "ts", "matrix")
  } else {
    class(y) <- "ts"
  }
  if("msts" %in% class(x))
    class(y) <- c("msts", class(y))
  attr <- attributes(x)
  attributes(y)$tsp <- attr$tsp
  attributes(y)$msts <- attr$msts
  return(y)
}

#' @export
print.scp <- function(x, ...) {
  NextMethod()
}

#' @export
summary.scp <- function(object, ...) {
  NextMethod()
}

#' @export
print.summary.scp <- function(x, ...) {
  NextMethod()
}

#' @export
print.cpforecast <- function(x, ...) {
  cat(paste(toupper(x$method), "\n\n"))
  if (!is.null(x$call)) {
    cat(paste("Call:\n"))
    for (i in 1:length(deparse(x$call))) {
      cat(paste("", deparse(x$call)[i]), "\n")
    }
    cat(paste("\n"))
  }

  cat(paste("", "cp_times =", x$cp_times,
            ifelse("mean" %in% names(x), "(the forward step included)", ""), "\n"))

  if ("model" %in% names(x)) {
    cat(paste("\nForecasts of the forward step:\n"))
    NextMethod()
  }
}

#' @export
summary.cpforecast <- function(object, ...) {
  class(object) <- c("summary.cpforecast", class(object))
  object
}

#' @export
print.summary.cpforecast <- function(x, ...) {
  NextMethod()
  cat("\nCross-validation error measures:\n")
  print(round(
    accuracy.default(x, measures = c(point_measures, interval_measures),
                     byhorizon = FALSE),
    digits = 3))
}
