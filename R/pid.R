#' Conformal PID control method
#'
#' Compute prediction intervals and other information by
#' applying the conformal PID (Proportional-Integral-Derivative) control method.
#'
#' The PID method combines three modules to make the final iteration:
#' \deqn{q_{t+h|t}=\underbrace{q_{t+h-1|t-1} + \eta(\mathrm{err}_{t|t-h}-\alpha)}_{\mathrm{P}}+\underbrace{r_t\left(\sum_{i=1}^t\left(\mathrm{err}_{i|i-h}-\alpha\right)\right)}_{\mathrm{I}}+\underbrace{\hat{s}_{t+h|t}}_{\mathrm{D}}}
#' for each individual forecast horizon \code{h}, respectively, where
#'   - Quantile tracking part (P) is \eqn{q_{t+h-1|t-1} + \eta(\mathrm{err}_{t|t-h}-\alpha)}, where \eqn{q_{1+h|1}} is set to 0 without a loss of generality, \eqn{\mathrm{err}_{t|t-h}=1} if \eqn{s_{t|t-h}>q_{t|t-h}}, and \eqn{\mathrm{err}_{t|t-h}=0} if \eqn{s_{t|t-h} \leq q_{t|t-h}}.
#'   - Error integration part (I) is \eqn{r_t\left(\sum_{i=1}^t\left(\mathrm{err}_{i|i-h}-\alpha\right)\right)}. Here we use a nonlinear saturation
#'   function \eqn{r_t(x)=K_{\mathrm{I}} \tan \left(x \log (t) /\left(t C_{\text {sat }}\right)\right)}, where we set \eqn{\tan (x)=\operatorname{sign}(x) \cdot \infty} for \eqn{x \notin[-\pi / 2, \pi / 2]}, and \eqn{C_{\text {sat }}, K_{\mathrm{I}}>0} are constants that we choose heuristically.
#'   - Scorecasting part (D) is \eqn{\hat{s}_{t+h|t}} is forecast generated
#'   by training a scorecaster based on nonconformity scores available at time \eqn{t}.
#'
#' @param object An object of class \code{"cvforecast"}. It must have an argument
#' \code{x} for original univariate time series, an argument \code{MEAN} for
#' point forecasts and \code{ERROR} for forecast errors on validation set.
#' See the results of a call to \code{\link{cvforecast}}.
#' @param alpha A numeric vector of significance levels to achieve a desired
#' coverage level \eqn{1-\alpha}. Defaults to the levels used in \code{object}.
#' @param symmetric If \code{TRUE}, symmetric nonconformity scores (i.e. \eqn{|e_{t+h|t}|})
#' are used. If \code{FALSE}, asymmetric nonconformity scores (i.e. \eqn{e_{t+h|t}})
#' are used, and then upper bounds and lower bounds are produced separately.
#' @param ncal Length of the burn-in period for training the scorecaster.
#' If \code{rolling = TRUE}, it is also used as the length of the trailing windows
#' for learning rate calculation and the windows for the calibration set.
#' If \code{rolling = FALSE}, it is used as initial period of calibration sets
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
#' @param lr A positive initial learning rate used for quantile tracking.
#' @param Tg The time that is set to achieve the target absolute coverage
#' guarantee before this. It must be greater than 1 when \code{Csat} is not
#' supplied. Defaults to the number of cross-validation periods in \code{object}.
#' @param delta A number in \eqn{(0, 1)}. The target absolute coverage guarantee
#' is set to \eqn{1-\alpha-\delta}.
#' @param Csat A positive constant ensuring that by time \code{Tg}, an absolute
#' guarantee is of at least \eqn{1-\alpha-\delta} coverage. Derived from \code{Tg}
#' and \code{delta} when not supplied.
#' @param KI A non-negative constant to place the integrator on the same scale as
#' the scores. Defaults to the largest absolute forecast error in \code{object}.
#' @param update If \code{TRUE}, \code{object} already holds the results of a
#' previous call and only the newly added time steps are computed; the
#' prediction intervals produced earlier are carried over unchanged. Set by
#' \code{\link{update.cpforecast}} and not normally set by hand.
#' @param ... Other arguments are passed to the \code{scorecastfun} function.
#'
#' @return A list of class \code{c("pid", "cpforecast", "cvforecast", "forecast")}
#' with the following components:
#' \item{x}{The original time series.}
#' \item{series}{The name of the series \code{x}.}
#' \item{xreg}{Exogenous predictor variables used, if applicable.}
#' \item{method}{A character string "pid".}
#' \item{cp_times}{An integer vector giving the number of conformal predictions
#' performed in cross-validation for each forecast horizon.}
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
#' \item{model}{A list containing information about the conformal prediction model.}
#' If \code{mean} is included in the \code{object}, the components \code{mean},
#' \code{lower}, and \code{upper} will also be returned, showing the information
#' about the forecasts generated using all available observations.
#'
#' @references Angelopoulos, A., Candes, E., and Tibshirani, R. J. (2024).
#' "Conformal PID control for time series prediction", \emph{Advances in Neural
#' Information Processing Systems}, \bold{36}, 23047--23074.
#' @examples
#' # Simulate time series from an AR(2) model
#' library(forecast)
#' set.seed(1)
#' series <- arima.sim(n = 200, list(ar = c(0.8, -0.5)), sd = sqrt(1))
#' # Cross-validation forecasting
#' far2 <- function(x, h, level) {
#'   Arima(x, order = c(2, 0, 0)) |>
#'     forecast(h = h, level)
#' }
#' fc <- cvforecast(series, forecastfun = far2, h = 3, level = 95,
#'                  forward = TRUE, initial = 1, window = 50)
#' # PID setup
#' Tg <- 200; delta <- 0.01
#' Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
#' KI <- 2
#' lr <- 0.1
#' # PID without scorecaster
#' pidfc_nsf <- pid(fc, symmetric = FALSE, ncal = 50, rolling = TRUE,
#'                  integrate = TRUE, scorecast = FALSE,
#'                  lr = lr, Tg = Tg, KI = KI, Csat = Csat)
#' print(pidfc_nsf)
#' summary(pidfc_nsf)
#' # PID with a Naive model for the scorecaster
#' naivefun <- function(x, h) {
#'   naive(x) |> forecast(h = h)
#' }
#' pidfc <- pid(fc, symmetric = FALSE, ncal = 50, rolling = TRUE,
#'              integrate = TRUE, scorecast = TRUE, scorecastfun = naivefun,
#'              lr = lr, Tg = Tg, KI = KI, Csat = Csat)
#'
#' @export
pid <- function(
  object,
  alpha = 1 - 0.01 * object$level,
  symmetric = FALSE,
  ncal = 10,
  rolling = FALSE,
  integrate = TRUE,
  scorecast = !symmetric,
  scorecastfun = NULL,
  lr = 0.1,
  Tg = NROW(object$ERROR),
  delta = 0.01,
  Csat = NULL,
  KI = max(abs(object$ERROR), na.rm = TRUE),
  update = FALSE,
  ...
) {
  # Check inputs
  if (any(alpha >= 1 | alpha <= 0)) {
    stop("alpha should be in (0, 1)")
  }
  if (ncal < 10) {
    stop("Length of calibration period should at least be 10")
  }
  if (scorecast && is.null(scorecastfun)) {
    stop("scorecastfun should not be NULL if scorecast is TRUE")
  }
  is_number <- function(x) {
    is.numeric(x) && length(x) == 1L && isTRUE(is.finite(x))
  }
  if (!is_number(lr) || lr <= 0) {
    stop("`lr` should be a positive number")
  }
  if (!is_number(KI) || KI < 0) {
    stop("`KI` should be a non-negative number")
  }

  if (is.null(Csat)) {
    if (!is_number(Tg) || Tg <= 1) {
      stop("`Tg` should be a number greater than 1")
    }
    if (!is_number(delta) || delta <= 0 || delta >= 1) {
      stop("`delta` should be in (0, 1)")
    }
    Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
  }
  if (!is_number(Csat) || Csat <= 0) {
    stop("`Csat` should be a positive number")
  }

  alpha <- sort(alpha, decreasing = TRUE)
  level <- 100 * (1 - alpha)
  nlev <- length(alpha)
  pf <- ts(
    as.matrix(object$MEAN),
    start = start(object$MEAN),
    frequency = frequency(object$MEAN)
  )
  errors <- ts(
    as.matrix(object$ERROR),
    start = start(object$ERROR),
    frequency = frequency(object$ERROR)
  )
  horizon <- ncol(pf)
  n <- nrow(pf)

  if (ncal + horizon - 1L > nrow(errors) - !object$forward) {
    stop(
      "`ncal` is too large: `ncal + h - 1` must not exceed ",
      nrow(errors) - !object$forward,
      " for h = ",
      horizon
    )
  }

  # Warm start: resume the recursion instead of replaying the whole history.
  # Only valid if the stored state exists and none of the settings changed.
  warm <- FALSE
  if (update && !is.null(object$model$state) && !is.null(object$model$t_last)) {
    settings <- list(
      alpha = alpha,
      symmetric = symmetric,
      ncal = ncal,
      rolling = rolling,
      integrate = integrate,
      scorecast = scorecast,
      lr = lr,
      Csat = Csat,
      KI = KI
    )
    unchanged <- vapply(
      names(settings),
      function(nm) {
        old <- object$model$args[[nm]]
        is.null(old) || isTRUE(all.equal(old, settings[[nm]]))
      },
      logical(1)
    )
    if (all(unchanged)) {
      warm <- TRUE
    } else {
      warning(
        "cannot warm-start `pid`: ",
        paste(names(settings)[!unchanged], collapse = ", "),
        " changed since the object was built; recomputing from scratch",
        call. = FALSE
      )
    }
  }

  namatrix <- `colnames<-`(
    ts(
      matrix(NA_real_, nrow = n, ncol = horizon),
      start = start(pf),
      frequency = frequency(pf)
    ),
    paste0("h=", seq(horizon))
  )
  nalist <- `names<-`(
    rep(list(namatrix), nlev),
    paste0(level, "%")
  )

  # Rows already covered by the stored state. Everything carried over from the
  # object is copied into the top of a freshly allocated full-size matrix, so
  # the time-series attributes come from `namatrix` and never have to be
  # rebuilt, and nothing has to be padded.
  prev <- integer(0)
  if (warm) {
    if (nrow(object$model$lr_update) > n) {
      stop("the stored state is longer than the current series")
    }
    prev <- seq_len(nrow(object$model$lr_update))
  }

  lower <- upper <- nalist
  lrmat <- namatrix
  if (warm) {
    for (i in seq(nlev)) {
      lower[[i]][seq_len(nrow(object$LOWER[[i]])), ] <- object$LOWER[[i]]
      upper[[i]][seq_len(nrow(object$UPPER[[i]])), ] <- object$UPPER[[i]]
    }
    lrmat[prev, ] <- object$model$lr_update
  }
  if (integrate) {
    integrator <- integrator_lower <- integrator_upper <- nalist
    if (warm) {
      for (i in seq(nlev)) {
        if (symmetric) {
          integrator[[i]][prev, ] <- object$model$integrator[[i]]
        } else {
          integrator_lower[[i]][prev, ] <- object$model$integrator$lower[[i]]
          integrator_upper[[i]][prev, ] <- object$model$integrator$upper[[i]]
        }
      }
    }
  }
  if (scorecast) {
    scorecaster <- scorecaster_lower <- scorecaster_upper <- namatrix
    if (warm) {
      if (symmetric) {
        scorecaster[prev, ] <- object$model$scorecaster
      } else {
        scorecaster_lower[prev, ] <- object$model$scorecaster$lower
        scorecaster_upper[prev, ] <- object$model$scorecaster$upper
      }
    }
  }

  out <- c(
    list(x = object$x, series = object$series),
    if ("xreg" %in% names(object)) list(xreg = object$xreg)
  )

  t_last <- nrow(errors) - !object$forward
  t_resume <- if (warm) object$model$t_last + 1L else NA_integer_

  cp_times <- integer(horizon)
  state <- vector("list", horizon)
  for (h in seq(horizon)) {
    idx_all <- seq(h, t_last, by = 1L)
    # number of times a conformal prediction has been made, over the whole
    # history -- not just over the steps replayed in this call
    cp_times[h] <- sum(idx_all >= ncal + h - 1)
    indx <- if (warm) idx_all[idx_all >= t_resume] else idx_all

    errt_h <- errt_lower_h <- errt_upper_h <-
      integ_h <- integ_lower_h <- integ_upper_h <-
        matrix(NA_real_, nrow = n, ncol = nlev)
    qts_h <- qts_lower_h <- qts_upper_h <-
      qs_h <- qs_lower_h <- qs_upper_h <-
        matrix(0, nrow = n, ncol = nlev)

    if (warm) {
      st <- object$model$state[[h]]
      if (symmetric) {
        qts_h[prev, ] <- st$qts
        qs_h[prev, ] <- st$qs
        errt_h[prev, ] <- st$errt
        qs_lower_h <- qs_upper_h <- qs_h
      } else {
        qts_lower_h[prev, ] <- st$qts_lower
        qts_upper_h[prev, ] <- st$qts_upper
        qs_lower_h[prev, ] <- st$qs_lower
        qs_upper_h[prev, ] <- st$qs_upper
        errt_lower_h[prev, ] <- st$errt_lower
        errt_upper_h[prev, ] <- st$errt_upper
      }
      # carry the integrator history over as well, so that writing it back
      # below does not overwrite the restored rows with NA
      if (integrate) {
        for (i in seq(nlev)) {
          if (symmetric) {
            integ_h[, i] <- integrator[[i]][, h]
          } else {
            integ_lower_h[, i] <- integrator_lower[[i]][, h]
            integ_upper_h[, i] <- integrator_upper[[i]][, h]
          }
        }
      }
    }

    for (t in indx) {
      t_burnin <- max(t - ncal + 1L, h)
      errors_subset <- subset(
        errors[, h],
        start = ifelse(!rolling, h, t_burnin),
        end = t
      )

      if (symmetric) {
        # Calculate errt
        errt_h[t, ] <- (abs(errors[t, h]) > qs_h[t, ])

        # Learning rate (same for the upper and lower bounds)
        lrmat[t, h] <- ifelse(
          length(errors_subset) <= 1,
          lr,
          lr *
            (max(abs(errors_subset), na.rm = TRUE) -
              min(abs(errors_subset), na.rm = TRUE))
        )

        # Update quantile tracking
        qts_h[t + h, ] <- qts_h[t + h - 1, ] +
          lrmat[t, h] * (errt_h[t, ] - alpha)

        # Update integrator
        if (integrate) {
          es <- errt_h[h:t, ] |> matrix(ncol = nlev)
          integrator_arg <- apply(es, 2, sum) - nrow(es) * alpha
          integ_h[t + h, ] <- sapply(
            seq(nlev),
            function(i) {
              ifelse(
                nrow(es) == 1,
                0,
                saturation_fn_log(integrator_arg[i], nrow(es), Csat, KI)
              )
            }
          )
        }

        # Update scorecaster
        do_scorecast <- (scorecast && t >= (ncal + h - 1))
        if (do_scorecast) {
          sc <- try(
            suppressWarnings(
              # h-step-ahead forecast
              scorecastfun(abs(errors_subset), h = h, ...)
            ),
            silent = TRUE
          )
          if (!is.element("try-error", class(sc))) {
            scorecaster[t + h, h] <- as.numeric(sc$mean[h])
          }
        }

        # Update the next quantile
        qs_h[t + h, ] <- qts_h[t + h, ] +
          ifelse(
            rep(integrate, nlev),
            integ_h[t + h, ],
            rep(0, nlev)
          ) +
          rep(
            ifelse(
              do_scorecast,
              ifelse(is.na(scorecaster[t + h, h]), 0, scorecaster[t + h, h]),
              0
            ),
            nlev
          )
        qs_lower_h[t + h, ] <- qs_upper_h[t + h, ] <- qs_h[t + h, ]
      } else {
        # Calculate errt
        errt_lower_h[t, ] <- (-errors[t, h]) > qs_lower_h[t, ]
        errt_upper_h[t, ] <- errors[t, h] > qs_upper_h[t, ]

        # Learning rate (same for the upper and lower bounds)
        lrmat[t, h] <- ifelse(
          length(errors_subset) <= 1,
          lr,
          lr *
            (max(errors_subset, na.rm = TRUE) -
              min(errors_subset, na.rm = TRUE))
        )

        # Update quantile tracking
        qts_lower_h[t + h, ] <- qts_lower_h[t + h - 1, ] +
          lrmat[t, h] * (errt_lower_h[t, ] - alpha / 2)
        qts_upper_h[t + h, ] <- qts_upper_h[t + h - 1, ] +
          lrmat[t, h] * (errt_upper_h[t, ] - alpha / 2)

        # Update integrator
        if (integrate) {
          el <- errt_lower_h[h:t, ] |> matrix(ncol = nlev)
          integrator_lower_arg <- apply(el, 2, sum) - nrow(el) * alpha / 2
          integ_lower_h[t + h, ] <- sapply(
            seq(nlev),
            function(i) {
              ifelse(
                nrow(el) == 1,
                0,
                saturation_fn_log(integrator_lower_arg[i], nrow(el), Csat, KI)
              )
            }
          )

          eu <- errt_upper_h[h:t, ] |> matrix(ncol = nlev)
          integrator_upper_arg <- apply(eu, 2, sum) - nrow(eu) * alpha / 2
          integ_upper_h[t + h, ] <- sapply(
            seq(nlev),
            function(i) {
              ifelse(
                nrow(eu) == 1,
                0,
                saturation_fn_log(integrator_upper_arg[i], nrow(eu), Csat, KI)
              )
            }
          )
        }

        # Update scorecaster
        do_scorecast <- (scorecast && t >= (ncal + h - 1))
        if (do_scorecast) {
          sc <- try(
            suppressWarnings(
              # h-step-ahead forecast
              scorecastfun(errors_subset, h = h, ...)
            ),
            silent = TRUE
          )
          if (!is.element("try-error", class(sc))) {
            scorecaster_lower[t + h, h] <- -as.numeric(sc$mean[h])
            scorecaster_upper[t + h, h] <- as.numeric(sc$mean[h])
          }
        }

        # Update the next quantile
        qs_lower_h[t + h, ] <- qts_lower_h[t + h, ] +
          ifelse(
            rep(integrate, nlev),
            integ_lower_h[t + h, ],
            rep(0, nlev)
          ) +
          rep(
            ifelse(
              do_scorecast,
              ifelse(
                is.na(scorecaster_lower[t + h, h]),
                0,
                scorecaster_lower[t + h, h]
              ),
              0
            ),
            nlev
          )
        qs_upper_h[t + h, ] <- qts_upper_h[t + h, ] +
          ifelse(
            rep(integrate, nlev),
            integ_upper_h[t + h, ],
            rep(0, nlev)
          ) +
          rep(
            ifelse(
              do_scorecast,
              ifelse(
                is.na(scorecaster_upper[t + h, h]),
                0,
                scorecaster_upper[t + h, h]
              ),
              0
            ),
            nlev
          )
      }

      # PIs
      if (t >= (ncal + h - 1)) {
        for (i in seq(nlev)) {
          lower[[i]][t + h, h] <- pf[t + h, h] - qs_lower_h[t + h, i]
          upper[[i]][t + h, h] <- pf[t + h, h] + qs_upper_h[t + h, i]
        }
      }
    }

    if (integrate) {
      for (i in seq(nlev)) {
        integrator[[i]][, h] <- integ_h[, i]
        integrator_lower[[i]][, h] <- integ_lower_h[, i]
        integrator_upper[[i]][, h] <- integ_upper_h[, i]
      }
    }
    state[[h]] <- if (symmetric) {
      list(qts = qts_h, qs = qs_h, errt = errt_h)
    } else {
      list(
        qts_lower = qts_lower_h,
        qts_upper = qts_upper_h,
        qs_lower = qs_lower_h,
        qs_upper = qs_upper_h,
        errt_lower = errt_lower_h,
        errt_upper = errt_upper_h
      )
    }
  }

  out$method <- paste("pid")
  out$cp_times <- cp_times
  out$MEAN <- object$MEAN
  out$ERROR <- object$ERROR
  out$LOWER <- lower
  out$UPPER <- upper
  out$level <- level
  out$call <- match.call()
  if ("mean" %in% names(object)) {
    out$mean <- object$mean
    out$lower <- extract_final(
      lower,
      nrow = n,
      ncol = horizon,
      bench = out$mean
    )
    out$upper <- extract_final(
      upper,
      nrow = n,
      ncol = horizon,
      bench = out$mean
    )
  }
  out$model$method <- out$method
  out$model$call <- match.call()
  # Every setting is recorded here and nowhere else; `call` is for display
  # only, and `update.cpforecast()` replays these values.
  out$model$args <- c(
    list(
      alpha = alpha,
      symmetric = symmetric,
      ncal = ncal,
      rolling = rolling,
      integrate = integrate,
      scorecast = scorecast,
      scorecastfun = scorecastfun,
      lr = lr,
      Tg = Tg,
      delta = delta,
      Csat = Csat,
      KI = KI
    ),
    list(...)
  )
  out$model$lr_update <- lrmat
  out$model$t_last <- t_last
  out$model$state <- state
  if (update) {
    out$model$cvforecast$call <- object$model$cvforecast$call
    out$model$cvforecast$args <- object$model$cvforecast$args
  } else {
    out$model$cvforecast$call <- object$call
    out$model$cvforecast$args <- object$args
  }
  if (symmetric) {
    if (integrate) {
      out$model$integrator <- integrator
    }
    if (scorecast) {
      out$model$scorecaster <- scorecaster
    }
  } else {
    if (integrate) {
      out$model$integrator <- list(
        lower = integrator_lower,
        upper = integrator_upper
      )
    }
    if (scorecast) {
      out$model$scorecaster <- list(
        lower = scorecaster_lower,
        upper = scorecaster_upper
      )
    }
  }

  return(structure(
    out,
    class = c("pid", "cpforecast", "cvforecast", "forecast")
  ))
}

saturation_fn_log <- function(x, t, Csat, KI) {
  if (KI == 0) {
    return(0)
  } else {
    tan_out <- mytan(x * log(t) / (Csat * (t)))
    out <- KI * tan_out
    return(out)
  }
}

mytan <- function(x) {
  if (x >= pi / 2) {
    return(Inf)
  } else if (x <= -pi / 2) {
    return(-Inf)
  } else {
    return(tan(x))
  }
}
