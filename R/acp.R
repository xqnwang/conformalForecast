#' Adaptive conformal prediction method
#'
#' Compute prediction intervals and other information by
#' applying the adaptive conformal prediction (ACP) method.
#'
#' The ACP method considers the online update:
#' \deqn{\alpha_{t+h|t}:=\alpha_{t+h-1|t-1}+\gamma(\alpha-\mathrm{err}_{t|t-h}),}
#' for each individual forecast horizon \code{h}, respectively,
#' where \eqn{\mathrm{err}_{t|t-h}=1} if \eqn{s_{t|t-h}>q_{t|t-h}}, and
#' \eqn{\mathrm{err}_{t|t-h}=0} if \eqn{s_{t|t-h} \leq q_{t|t-h}}.
#'
#' @aliases print.acp summary.acp print.summary.acp
#'
#' @param object An object of class \code{"cvforecast"}. It must have an argument
#' \code{x} for original univariate time series, an argument \code{MEAN} for
#' point forecasts and \code{ERROR} for forecast errors on validation set.
#' See the results of a call to \code{\link{cvforecast}}.
#' @param alpha A numeric vector of significance levels to achieve a desired
#' coverage level \eqn{1-\alpha}.
#' @param gamma The step size parameter \eqn{\gamma>0} for \eqn{\alpha} updating.
#' @param symmetric If \code{TRUE}, symmetric nonconformity scores (i.e. \eqn{|e_{t+h|t}|})
#' are used. If \code{FALSE}, asymmetric nonconformity scores (i.e. \eqn{e_{t+h|t}})
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
#' @param na.rm If \code{TRUE}, corresponding entries in sample values are removed
#' if it is \code{NA} when calculating sample quantile.
#' @param ... Other arguments are passed to the
#' \code{\link[ggdist]{weighted_quantile}} function for quantile computation.
#'
#' @return A list of class \code{c("acp", "cpforecast", "forecast")}
#' with the following components:
#' \item{x}{The original time series.}
#' \item{series}{The name of the series \code{x}.}
#' \item{method}{A character string "acp".}
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
#' If \code{mean} is included in the \code{object}, the components \code{mean},
#' \code{lower}, and \code{upper} will also be returned, showing the information
#' about the forecasts generated using all available observations.
#'
#' @references Gibbs, I., and Candes, E. (2021). "Adaptive conformal inference under
#' distribution shift", \emph{Advances in Neural Information Processing Systems},
#' \bold{34}, 1660--1672.
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
#' # ACP with asymmetric nonconformity scores and rolling calibration sets
#' acpfc <- acp(fc, symmetric = FALSE, gamma = 0.005, ncal = 100, rolling = TRUE)
#' print(acpfc)
#' summary(acpfc)
#'
#' @importFrom ggdist weighted_quantile
#' @export
acp <- function(object, alpha = 1 - 0.01 * object$level, gamma = 0.005,
                symmetric = FALSE, ncal = 10, rolling = FALSE,
                quantiletype = 1, na.rm = TRUE, ...) {
  # Check inputs
  if (any(alpha >= 1 | alpha <= 0))
    stop("alpha should be in (0, 1)")
  if (gamma < 0)
    stop("the step size parameter gamma should be positive")
  if (ncal < 10)
    stop("length of calibration period should at least be 10")
  if (!quantiletype %in% 1:9)
    stop("quantiletype is invalid. It must be in 1:9.")

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
  lower <- upper <- alphat <- alphat_lower <- alphat_upper <-
    `names<-` (rep(list(namatrix), length(alpha)),
               paste0(level, "%"))

  out <- list(
    x = object$x,
    series = object$series
  )

  for (h in seq(horizon)) {
    indx <- seq(ncal+h-1, nrow(errors), by = 1L)

    alphat_h <- alphat_lower_h <- alphat_upper_h <-
      errt_h <- errt_lower_h <- errt_upper_h <-
      q_lo_h <- q_up_h <-
      matrix(NA_real_, nrow = n, ncol = length(alpha))

    for (t in indx) {
      errors_subset <- subset(
        errors[, h],
        start = ifelse(!rolling, 1, t - ncal + 1L),
        end = t)

      if (symmetric) {
        if (t == indx[1])
          alphat_h[t+h, ] <- alpha

        # Compute sample quantiles
        q_lo_h[t+h, ] <- q_up_h[t+h, ] <- ggdist::weighted_quantile(
          x = abs(c(errors_subset, Inf)),
          probs = 1 - alphat_h[t+h, ],
          type = quantiletype,
          na.rm = na.rm,
          ...)

        # Compute errt
        tryCatch(
          {
            errt_h[t+1, ] <- abs(errors[t+1, h]) > q_lo_h[t+1, ]
            outl <- which(alphat_h[t+1, ] >= 1)
            outs <- which(alphat_h[t+1, ] <= 0)
            errt_h[t+1, outl] <- TRUE
            errt_h[t+1, outs] <- FALSE
          },
          error = function(e) {
            errt_h[t+1, ] <- NA_real_
          }
        )

        if (t < tail(indx, 1)) {
          if (any(is.na(errt_h[t+1, ]))) {
            # Keep alpha unchanged
            alphat_h[t+h+1, ] <- alphat_h[t+h, ]
          } else {
            # Update alpha
            alphat_h[t+h+1, ] <- alphat_h[t+h, ] + gamma*(alpha - errt_h[t+1, ])
          }
        }
      } else {
        if (t == indx[1])
          alphat_lower_h[t+h, ] <- alphat_upper_h[t+h, ] <- alpha/2

        # Compute sample quantiles
        q_lo_h[t+h, ] <- ggdist::weighted_quantile(
          x = -c(errors_subset, Inf),
          probs = 1 - alphat_lower_h[t+h, ],
          type = quantiletype,
          na.rm = na.rm,
          ...)
        q_up_h[t+h, ] <- ggdist::weighted_quantile(
          x = c(errors_subset, Inf),
          probs = 1 - alphat_upper_h[t+h, ],
          type = quantiletype,
          na.rm = na.rm,
          ...)

        # Compute errt
        tryCatch(
          {
            errt_lower_h[t+1, ] <- (-errors[t+1, h]) > q_lo_h[t+1, ]
            errt_lower_h[t+1, which(alphat_lower_h[t+1, ] >= 1)] <- TRUE
            errt_lower_h[t+1, which(alphat_lower_h[t+1, ] <= 0)] <- FALSE

            errt_upper_h[t+1, ] <- (errors[t+1, h]) > q_up_h[t+1, ]
            errt_upper_h[t+1, which(alphat_upper_h[t+1, ] >= 1)] <- TRUE
            errt_upper_h[t+1, which(alphat_upper_h[t+1, ] <= 0)] <- FALSE
          },
          error = function(e) {
            errt_lower_h[t+1, ] <- NA_real_
            errt_upper_h[t+1, ] <- NA_real_
          }
        )

        if (t < tail(indx, 1)) {
          if (any(is.na(errt_lower_h[t+1, ])) || any(is.na(errt_upper_h[t+1, ]))) {
            # Keep alpha unchanged
            alphat_lower_h[t+h+1, ] <- alphat_lower_h[t+h, ]
            alphat_upper_h[t+h+1, ] <- alphat_upper_h[t+h, ]
          } else {
            # Update alpha
            alphat_lower_h[t+h+1, ] <- alphat_lower_h[t+h, ] + gamma*(alpha/2 - errt_lower_h[t+1, ])
            alphat_upper_h[t+h+1, ] <- alphat_upper_h[t+h, ] + gamma*(alpha/2 - errt_upper_h[t+1, ])
          }
        }
      }
      for (i in seq(length(alpha))) {
        lower[[i]][t+h, h] <- pf[t+h, h] - q_lo_h[t+h, i]
        upper[[i]][t+h, h] <- pf[t+h, h] + q_up_h[t+h, i]
      }
    }
    for (i in seq(length(alpha))) {
      alphat[[i]][, h] <- alphat_h[, i]
      alphat_lower[[i]][, h] <- alphat_lower_h[, i]
      alphat_upper[[i]][, h] <- alphat_upper_h[, i]
    }
  }

  out$method <- paste("acp")
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
  out$model$gamma <- gamma
  out$model$symmetric <- symmetric
  if (symmetric) {
    out$model$alpha_update <- alphat
  } else {
    out$model$alpha_update <- list(lower = alphat_lower, upper = alphat_upper)
  }

  return(structure(out, class = c("acp", "cpforecast", "forecast")))
}

#' @export
print.acp <- function(x, ...) {
  NextMethod()
}

#' @export
summary.acp <- function(object, ...) {
  NextMethod()
}

#' @export
print.summary.acp <- function(x, ...) {
  NextMethod()
}
