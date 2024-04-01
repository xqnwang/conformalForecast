#' Conformal PID control method
#'
#' Compute prediction intervals and other information obtained by
#' applying the conformal PID control method.
#'
#' The PID method combines three modules to make the final iteration:
#' \deqn{q_{t+h|t}=\underbrace{p_{t+h-1|t-1} + \eta(\mathrm{err}_{t|t-h}-\alpha)}_{\mathrm{P}}+\underbrace{r_t\left(\sum_{i=1}^t\left(\mathrm{err}_{i|i-h}-\alpha\right)\right)}_{\mathrm{I}}+\underbrace{\hat{s}_{t+h|t}}_{\mathrm{D}}}
#' for each individual forecast horizon \code{h}, respectively, where
#'   - Quantile tracking part (P) is \eqn{p_{t+h-1|t-1} + \eta(\mathrm{err}_{t|t-h}-\alpha)}, where \eqn{p_{1+h|1}} is set to 0 without a loss of generality.
#'   - Error integration part (I) is \eqn{r_t\left(\sum_{i=1}^t\left(\mathrm{err}_{i|i-h}-\alpha\right)\right)}. Here we use a nonlinear saturation
#'   function \eqn{r_t(x)=K_{\mathrm{I}} \tan \left(x \log (t) /\left(t C_{\text {sat }}\right)\right)}, where we set \eqn{\tan (x)=\operatorname{sign}(x) \cdot \infty} for \eqn{x \notin[-\pi / 2, \pi / 2]}, and \eqn{C_{\text {sat }}, K_{\mathrm{I}}>0} are constants that we choose heuristically.
#'   - Scorecasting part (D) is \eqn{\hat{s}_{t+h|t}} is forecast generated
#'   by training a scorecaster based on conformity scores available at time \eqn{t}.
#'
#' @aliases print.pid summary.pid print.summary.pid
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
#' @param ncal Length of the burn-in period for training the scorecaster.
#' If \code{rolling = FALSE}, it also denotes the length of the trailing windows
#' for learning rate calculation and the windows for the calibration set.
#' If \code{rolling = FALSE}, it denotes the initial period of calibration sets
#' and trailing windows for learning rate calculation.
#' @param rolling If \code{TRUE}, a rolling window strategy will be adopted to
#' form the trailing window for learning rate calculation and the calibration set
#' for scorecaster if applicable. Otherwise, expanding window strategy will be used.
#' @param integrate If \code{TRUE}, error integration will be included in the
#' update process.
#' @param scorecast If \code{TRUE}, scorecasting will be included in the update
#' process, and \code{scorecastfun} should be given.
#' @param scorecastfun A scorecaster function to return an object of class
#' \code{forecast}. Its first argument must be a univariate time series, and
#' it must have an argument \code{h} for the forecast horizon.
#' @param lr Initial learning rate used for quantile tracking.
#' @param Tg The time that is set to achieve the target absolute coverage
#' guarantee before this.
#' @param delta The target absolute coverage guarantee is set to \eqn{1-\alpha-\delta}.
#' @param Csat A positive constant ensuring that by time \code{Tg}, an absolute
#' guarantee is of at least \eqn{1-\alpha-\delta} coverage.
#' @param KI A positive constant to place the integrator on the same scale as the scores.
#' @param ... Other arguments are passed to the \code{scorecastfun} function.
#'
#' @return A list of class \code{c("pid", "cpforecast", "forecast")}
#' with the following components:
#' \item{x}{The original time series.}
#' \item{series}{The name of the series \code{x}.}
#' \item{method}{A character string "pid".}
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
#' @references Angelopoulos, A., Candes, E., and Tibshirani, R. J. (2024).
#' "Conformal PID control for time series prediction", \emph{Advances in Neural
#' Information Processing Systems}, \bold{36}, 23047--23074.
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
#' # PID setup
#' Tg <- 1000; delta <- 0.01
#' Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
#' KI <- 2
#' lr <- 0.1
#'
#' # PID without scorecaster
#' pidfc_nsf <- pid(fc, symmetric = FALSE, ncal = 100, rolling = TRUE,
#'                  integrate = TRUE, scorecast = FALSE,
#'                  lr = lr, KI = KI, Csat = Csat)
#' print(pidfc_nsf)
#' summary(pidfc_nsf)
#'
#' # PID with a Theta model for the scorecaster
#' naivefun <- function(x, h) {
#'   naive(x) |> forecast(h = h)
#' }
#' pidfc <- pid(fc, symmetric = FALSE, ncal = 100, rolling = TRUE,
#'              integrate = TRUE, scorecast = TRUE, scorecastfun = naivefun,
#'              lr = lr, KI = KI, Csat = Csat)
#'
#' @export
pid <- function(object, alpha = 1 - 0.01 * object$level,
                symmetric = FALSE, ncal = 10, rolling = FALSE,
                integrate = TRUE, scorecast = !symmetric, scorecastfun = NULL,
                lr = 0.1, Tg = NULL, delta = NULL,
                Csat = 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg)),
                KI = abs(object$errors) |> max(na.rm = TRUE), ...) {
  # Check inputs
  if (any(alpha >= 1 | alpha <= 0))
    stop("alpha should be in (0, 1)")
  if (ncal < 10)
    stop("Length of calibration period should at least be 10")
  if (scorecast && is.null(scorecastfun))
    stop("scorecastfun should not be NULL if scorecast is TRUE")
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
    integrator <- integrator_lower <- integrator_upper <- nalist
  if (scorecast)
    scorecaster <- scorecaster_lower <- scorecaster_upper <- namatrix

  out <- list(
    x = object$x,
    series = object$series
  )

  for (h in seq(horizon)) {
    indx <- seq(h, nrow(errors), by = 1L)

    errt_h <- errt_lower_h <- errt_upper_h <-
      integ_h <- integ_lower_h <- integ_upper_h <-
      q_lo_h <- q_up_h <-
      matrix(NA_real_, nrow = n, ncol = length(alpha))
    qts_h <- qts_lower_h <- qts_upper_h <-
      qs_h <- qs_lower_h <- qs_upper_h <-
      matrix(0, nrow = n, ncol = length(alpha))

    for (t in indx) {
      t_burnin <- max(t - ncal + 1L, h)
      errors_subset <- subset(
        errors[, h],
        start = ifelse(!rolling, h, t_burnin),
        end = t)

      if (symmetric) {
        # Calculate errt
        errt_h[t, ] <- (abs(errors[t, h]) > qs_h[t, ])

        # Learning rate (same for the upper and lower bounds)
        lrmat[t, h] <- ifelse(
          length(errors_subset) <= 1,
          lr,
          lr*(max(abs(errors_subset), na.rm = TRUE) - min(abs(errors_subset), na.rm = TRUE)))

        # Update quantile tracking
        qts_h[t+h, ] <- qts_h[t+h-1, ] + lrmat[t, h] * (errt_h[t, ] - alpha)

        # Update integrator
        if (integrate) {
          es <- errt_h[h:t, ] |> matrix(ncol = length(alpha))
          integrator_arg <- apply(es, 2, sum) - nrow(es)*alpha
          integ_h[t+h, ] <- sapply(
            1:length(alpha),
            function(i) ifelse(
              nrow(es) == 1,
              0,
              saturation_fn_log(integrator_arg[i], nrow(es), Csat, KI)))
        }

        # Update scorecaster
        do_scorecast <- (scorecast && t >= (ncal+h-1))
        if (do_scorecast) {
          sc <- try(suppressWarnings(
            # h-step-ahead forecast
            scorecastfun(abs(errors_subset), h = h, ...)
          ), silent = TRUE)
          if (!is.element("try-error", class(sc))) {
            scorecaster[t+h, h] <- as.numeric(sc$mean[h])
          }
        }

        # Update the next quantile
        qs_h[t+h, ] <- qts_h[t+h, ] +
          ifelse(rep(integrate, length(alpha)), integ_h[t+h, ], rep(0, length(alpha))) +
          rep(ifelse(do_scorecast, scorecaster[t+h, h], 0), length(alpha))
        qs_lower_h[t+h, ] <- qs_upper_h[t+h, ] <- qs_h[t+h, ]
      } else {
        # Calculate errt
        errt_lower_h[t, ] <- (-errors[t, h]) > qs_lower_h[t, ]
        errt_upper_h[t, ] <- errors[t, h] > qs_upper_h[t, ]

        # Learning rate (same for the upper and lower bounds)
        lrmat[t, h] <- ifelse(
          length(errors_subset) <= 1,
          lr,
          lr*(max(errors_subset, na.rm = TRUE) - min(errors_subset, na.rm = TRUE)))

        # Update quantile tracking
        qts_lower_h[t+h, ] <- qts_lower_h[t+h-1, ] + lrmat[t, h] * (errt_lower_h[t, ] - alpha/2)
        qts_upper_h[t+h, ] <- qts_upper_h[t+h-1, ] + lrmat[t, h] * (errt_upper_h[t, ] - alpha/2)

        # Update integrator
        if (integrate) {
          el <- errt_lower_h[h:t, ] |> matrix(ncol = length(alpha))
          integrator_lower_arg <- apply(el, 2, sum) - nrow(el)*alpha/2
          integ_lower_h[t+h, ] <- sapply(
            1:length(alpha),
            function(i) ifelse(
              nrow(el) == 1,
              0,
              saturation_fn_log(integrator_lower_arg[i], nrow(el), Csat, KI)))

          eu <- errt_upper_h[h:t, ] |> matrix(ncol = length(alpha))
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
          sc <- try(suppressWarnings(
            # h-step-ahead forecast
            scorecastfun(errors_subset, h = h, ...)
          ), silent = TRUE)
          if (!is.element("try-error", class(sc))) {
            scorecaster_lower[t+h, h] <- - as.numeric(sc$mean[h])
            scorecaster_upper[t+h, h] <- as.numeric(sc$mean[h])
          }
        }

        # Update the next quantile
        qs_lower_h[t+h, ] <- qts_lower_h[t+h, ] +
          ifelse(rep(integrate, length(alpha)), integ_lower_h[t+h, ], rep(0, length(alpha))) +
          rep(ifelse(do_scorecast, scorecaster_lower[t+h, h], 0), length(alpha))
        qs_upper_h[t+h, ] <- qts_upper_h[t+h, ] +
          ifelse(rep(integrate, length(alpha)), integ_upper_h[t+h, ], rep(0, length(alpha))) +
          rep(ifelse(do_scorecast, scorecaster_upper[t+h, h], 0), length(alpha))
      }

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
        integrator[[i]][, h] <- integ_h[, i]
        integrator_lower[[i]][, h] <- integ_lower_h[, i]
        integrator_upper[[i]][, h] <- integ_upper_h[, i]
      }
    }
  }

  out$method <- paste("pid")
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
  out$model$integrate <- integrate
  out$model$scorecast <- scorecast
  out$model$lr <- lr
  out$model$Csat <- Csat
  out$model$KI <- KI
  out$model$lr_update <- lrmat
  if (symmetric) {
    if (integrate)
      out$model$integrator <- integrator
    if (scorecast)
      out$model$scorecaster <- scorecaster
  } else {
    if (integrate)
      out$model$integrator <- list(lower = integrator_lower, upper = integrator_upper)
    if (scorecast)
      out$model$scorecaster <- list(lower = scorecaster_lower, upper = scorecaster_upper)
  }

  return(structure(out, class = c("pid", "cpforecast", "forecast")))
}

saturation_fn_log <- function(x, t, Csat, KI) {
  if (KI == 0) {
    return(0)
  } else {
    tan_out <- mytan(x * log(t)/(Csat * (t)))
    out <- KI * tan_out
    return(out)
  }
}

mytan <- function(x){
  if (x >= pi/2) {
    return(Inf)
  } else if (x <= - pi/2) {
    return(-Inf)
  } else {
    return(tan(x))
  }
}

#' @export
print.pid <- function(x, ...) {
  NextMethod()
}

#' @export
summary.pid <- function(x, ...) {
  NextMethod()
}

#' @export
print.summary.pid <- function(x, ...) {
  NextMethod()
}
