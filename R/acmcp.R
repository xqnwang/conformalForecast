#' Autocorrelated multistep-ahead conformal prediction method
#'
#' Compute prediction intervals and other information by applying the Autocorrelated
#' Multistep-ahead Conformal Prediction (AcMCP) method. The method can only
#' deal with asymmetric nonconformity scores, i.e., forecast errors.
#'
#' Similar to the PID method, the AcMCP method also integrates three modules (P, I, and D) to
#' form the final iteration. However, instead of performing conformal prediction
#' for each individual forecast horizon \code{h} separately, AcMCP employs a combination
#' of an MA\eqn{(h-1)} model and a linear regression model of \eqn{e_{t+h|t}} on
#' \eqn{e_{t+h-1|t},\dots,e_{t+1|t}} as the scorecaster. This allows the AcMCP method
#' to capture the relationship between the \eqn{h}-step ahead forecast error and
#' past errors.
#'
#' @aliases print.acmcp summary.acmcp print.summary.acmcp
#'
#' @param object An object of class \code{"cvforecast"}. It must have an argument
#' \code{x} for original univariate time series, an argument \code{MEAN} for
#' point forecasts and \code{ERROR} for forecast errors on validation set.
#' See the results of a call to \code{\link{cvforecast}}.
#' @param alpha A numeric vector of significance levels to achieve a desired
#' coverage level \eqn{1-\alpha}.
#' @param ncal Length of the burn-in period for training the scorecaster.
#' If \code{rolling = TRUE}, it is also used as the length of the trailing windows
#' for learning rate calculation and the windows for the calibration set.
#' If \code{rolling = FALSE}, it is used as the initial period of calibration sets
#' and trailing windows for learning rate calculation.
#' @param rolling If \code{TRUE}, a rolling window strategy will be adopted to
#' form the trailing window for learning rate calculation and the calibration set
#' for scorecaster if applicable. Otherwise, expanding window strategy will be used.
#' @param integrate If \code{TRUE}, error integration will be included in the
#' update process.
#' @param scorecast If \code{TRUE}, scorecasting will be included in the update
#' process.
#' @param lr Initial learning rate used for quantile tracking.
#' @param Tg The time that is set to achieve the target absolute coverage
#' guarantee before this.
#' @param delta The target absolute coverage guarantee is set to \eqn{1-\alpha-\delta}.
#' @param Csat A positive constant ensuring that by time \code{Tg}, an absolute
#' guarantee is of at least \eqn{1-\alpha-\delta} coverage.
#' @param KI A positive constant to place the integrator on the same scale as the scores.
#' @param ... Other arguments are passed to the function.
#'
#' @return A list of class \code{c("acmcp", "cpforecast", "forecast")}
#' with the following components:
#' \item{x}{The original time series.}
#' \item{series}{The name of the series \code{x}.}
#' \item{method}{A character string "acmcp".}
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
#' about the test set forecasts generated using all available observations.
#'
#' @references Wang, X., and Hyndman, R. J. (2024). "Online conformal inference
#' for multi-step time series forecasting", arXiv preprint arXiv:2410.13115.
#' @seealso \code{\link{pid}}
#' @examples
#' # Simulate time series from an AR(2) model
#' library(forecast)
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
#' # AcMCP setup
#' Tg <- 200; delta <- 0.01
#' Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
#' KI <- 2
#' lr <- 0.1
#'
#' # AcMCP with integrator and scorecaster
#' acmcpfc <- acmcp(fc, ncal = 50, rolling = TRUE,
#'              integrate = TRUE, scorecast = TRUE,
#'              lr = lr, KI = KI, Csat = Csat)
#' print(acmcpfc)
#' summary(acmcpfc)
#'
#' @importFrom stats lm
#' @importFrom forecast meanf Arima forecast
#' @export
acmcp <- function(object, alpha = 1 - 0.01 * object$level,
                ncal = 10, rolling = FALSE,
                integrate = TRUE, scorecast = TRUE,
                lr = 0.1, Tg = NULL, delta = NULL,
                Csat = 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg)),
                KI = abs(object$errors) |> max(na.rm = TRUE), ...) {
  # Check inputs
  if (any(alpha >= 1 | alpha <= 0))
    stop("alpha should be in (0, 1)")
  if (ncal < 10)
    stop("Length of calibration period should at least be 10")
  if (!is.null(Tg) && !is.null(delta)) {
    if (!is.null(Csat))
      warning("Csat is replaced by calculation using Tg and delta")
    Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
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

  namatrix <- `colnames<-` (
    ts(matrix(NA_real_, nrow = n, ncol = horizon),
       start = start(pf),
       frequency = frequency(pf)),
    paste0("h=", seq(horizon)))
  nalist <- `names<-` (
    rep(list(namatrix), length(alpha)),
    paste0(level, "%"))

  lower <- upper <- nalist
  lrmat <- namatrix
  if (integrate)
    integrator_lower <- integrator_upper <- nalist
  if (scorecast)
    scorecaster_lower <- scorecaster_upper <- namatrix

  out <- list(
    x = object$x,
    series = object$series
  )

  for (h in seq(horizon)) {
    indx <- seq(h, nrow(errors)-!object$forward, by = 1L)

    errt_lower_h <- errt_upper_h <-
      integ_lower_h <- integ_upper_h <-
      q_lo_h <- q_up_h <-
      matrix(NA_real_, nrow = n, ncol = length(alpha))
    qts_lower_h <- qts_upper_h <-
      qs_lower_h <- qs_upper_h <- matrix(0, nrow = n, ncol = length(alpha))

    for (t in indx) {
      t_burnin <- max(t - ncal + 1L, h)
      errors_subset <- subset(
        errors[, h],
        start = ifelse(!rolling, h, t_burnin),
        end = t)

      # Calculate errt
      errt_lower_h[t, ] <- (-errors[t, h]) > qs_lower_h[t, ]
      errt_upper_h[t, ] <- errors[t, h] > qs_upper_h[t, ]

      # Learning rate (same for the upper and lower bounds)
      lrmat[t, h] <- ifelse(
        length(errors_subset) <= 1,
        lr,
        lr*(max(errors_subset) - min(errors_subset)))

      # Update quantile tracking
      qts_lower_h[t+h, ] <- qts_lower_h[t+h-1, ] + lrmat[t, h] * (errt_lower_h[t, ] - alpha/2)
      qts_upper_h[t+h, ] <- qts_upper_h[t+h-1, ] + lrmat[t, h] * (errt_upper_h[t, ] - alpha/2)

      # Update integrator
      if (integrate) {
        el <- errt_lower_h[h:t, ] |>
          matrix(ncol = length(alpha))
        integrator_lower_arg <- apply(el, 2, sum) - nrow(el)*alpha/2
        integ_lower_h[t+h, ] <- sapply(
          1:length(alpha),
          function(i) ifelse(
            nrow(el) == 1,
            0,
            saturation_fn_log(integrator_lower_arg[i], nrow(el), Csat, KI)))

        eu <- errt_upper_h[h:t, ] |>
          matrix(ncol = length(alpha))
        integrator_upper_arg <- apply(eu, 2, sum) - nrow(eu)*alpha/2
        integ_upper_h[t+h, ] <- sapply(
          1:length(alpha),
          function(i) ifelse(
            nrow(eu) == 1,
            0,
            saturation_fn_log(integrator_upper_arg[i], nrow(eu), Csat, KI)))
      }

      # Update scorecaster
      do_scorecast <- (scorecast && t >= (ncal+h-1))
      if (do_scorecast) {
        if (h == 1) {
          model <- forecast::meanf(errors_subset, h = h)
          scorecaster_lower[t+h, h] <- -as.numeric(model$mean)
          scorecaster_upper[t+h, h] <- as.numeric(model$mean)
        } else {
          model_MA <- forecast::Arima(errors_subset, order = c(0, 0, h-1)) |>
            forecast::forecast(h = h)
          model_LR <- lm(
            as.formula(paste0("V", h, " ~ .")),
            data = setNames(
              as.data.frame(sapply(1:h, function(j) {
                subset(
                  errors[, j],
                  start = ifelse(!rolling, j, t-ncal+1-h+j),
                  end = t-h+j)
              })),
              paste0("V", 1:h))
          ) |>
            forecast::forecast(
              newdata = setNames(
                data.frame(sapply(1:(h-1),
                                  function(j) scorecaster_upper[t-h+1+j, j]) |> matrix(nrow = 1)),
                paste0("V", 1:(h-1))))
          scorecaster_lower[t+h, h] <- -(as.numeric(model_MA$mean[h]) + as.numeric(model_LR$mean))/2
          scorecaster_upper[t+h, h] <- (as.numeric(model_MA$mean[h]) + as.numeric(model_LR$mean))/2
        }
      }

      # Update the next quantile
      qs_lower_h[t+h, ] <- qts_lower_h[t+h, ] +
        ifelse(rep(integrate, length(alpha)), integ_lower_h[t+h, ], rep(0, length(alpha))) +
        rep(ifelse(do_scorecast,
                   ifelse(is.na(scorecaster_lower[t+h, h]), 0, scorecaster_lower[t+h, h]),
                   0),
            length(alpha))
      qs_upper_h[t+h, ] <- qts_upper_h[t+h, ] +
        ifelse(rep(integrate, length(alpha)), integ_upper_h[t+h, ], rep(0, length(alpha))) +
        rep(ifelse(do_scorecast,
                   ifelse(is.na(scorecaster_upper[t+h, h]), 0, scorecaster_upper[t+h, h]),
                   0),
            length(alpha))

      # PIs
      if (t >= (ncal+h-1)) {
        for (i in seq(length(alpha))) {
          lower[[i]][t+h, h] <- pf[t+h, h] - qs_lower_h[t+h, i]
          upper[[i]][t+h, h] <- pf[t+h, h] + qs_upper_h[t+h, i]
        }
      }
    }
    if (integrate) {
      for (i in seq(length(alpha))) {
        integrator_lower[[i]][, h] <- integ_lower_h[, i]
        integrator_upper[[i]][, h] <- integ_upper_h[, i]
      }
    }
  }

  out$method <- paste("acmcp")
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
  out$model$integrate <- integrate
  out$model$scorecast <- scorecast
  out$model$lr <- lr
  out$model$Csat <- Csat
  out$model$KI <- KI
  out$model$lr_update <- lrmat
  if (integrate)
    out$model$integrator <- list(lower = integrator_lower, upper = integrator_upper)
  if (scorecast)
    out$model$scorecaster <- list(lower = scorecaster_lower, upper = scorecaster_upper)

  return(structure(out, class = c("acmcp", "cpforecast", "forecast")))
}


#' @export
print.acmcp <- function(x, ...) {
  NextMethod()
}

#' @export
summary.acmcp <- function(object, ...) {
  NextMethod()
}

#' @export
print.summary.acmcp <- function(x, ...) {
  NextMethod()
}
