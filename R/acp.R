#' Adaptive conformal prediction method
#'
#' Compute prediction intervals and other information by applying the adaptive
#' conformal prediction (ACP) method.
#'
#' The ACP method considers the online update:
#' \deqn{\alpha_{t+h|t}:=\alpha_{t+h-1|t-1}+\gamma(\alpha-\mathrm{err}_{t|t-h}),}
#' for each individual forecast horizon \code{h}, respectively, where
#' \eqn{\mathrm{err}_{t|t-h}=1} if \eqn{s_{t|t-h}>q_{t|t-h}}, and
#' \eqn{\mathrm{err}_{t|t-h}=0} if \eqn{s_{t|t-h} \leq q_{t|t-h}}.
#'
#' @inheritParams scp
#' @param gamma The step size parameter \eqn{\gamma>0} for the \eqn{\alpha}
#'   updating. Defaults to \code{0.005}.
#' @param na.rm If \code{TRUE}, corresponding entries in the sample values are
#'   removed if they are \code{NA} when calculating the sample quantile.
#'   Defaults to \code{TRUE}.
#' @param ... Other arguments are passed to the
#'   \code{\link[ggdist]{weighted_quantile}} function for the quantile
#'   computation.
#'
#' @return A list of class
#'   \code{c("acp", "cpforecast", "cvforecast", "forecast")} with the
#'   following components:
#'   \item{x}{The original time series.}
#'   \item{series}{The name of the series \code{x}.}
#'   \item{xreg}{Exogenous predictor variables used, if applicable.}
#'   \item{method}{A character string "acp".}
#'   \item{cp_times}{An integer vector giving the number of conformal
#'     predictions performed in cross-validation for each forecast horizon.}
#'   \item{MEAN}{Point forecasts as a multivariate time series, where the
#'     \eqn{h}th column holds the point forecasts for forecast horizon
#'     \eqn{h}. The time index corresponds to the period for which the
#'     forecast is produced.}
#'   \item{ERROR}{Forecast errors given by
#'     \eqn{e_{t+h|t} = y_{t+h}-\hat{y}_{t+h|t}}{e[t+h] = y[t+h]-f[t+h]}.}
#'   \item{LOWER}{A list containing lower bounds for prediction intervals for
#'     each \code{level}. Each element within the list will be a multivariate
#'     time series with the same dimensional characteristics as \code{MEAN}.}
#'   \item{UPPER}{A list containing upper bounds for prediction intervals for
#'     each \code{level}. Each element within the list will be a multivariate
#'     time series with the same dimensional characteristics as \code{MEAN}.}
#'   \item{level}{The confidence values associated with the prediction
#'     intervals.}
#'   \item{call}{The matched call.}
#'   \item{model}{A list containing information about the conformal prediction
#'     model: the resolved arguments in \code{model$args}, the call and the
#'     arguments of the underlying cross-validation in
#'     \code{model$cvforecast}, and the sequence of updated significance
#'     levels in \code{model$alpha_update}. The latter holds a single element
#'     \code{alpha} when \code{symmetric = TRUE}, and the two elements
#'     \code{lower} and \code{upper} when \code{symmetric = FALSE}; each is a
#'     list with one multivariate time series per confidence level, laid out
#'     like \code{MEAN}.}
#'   If \code{mean} is included in the \code{object}, the components
#'   \code{mean}, \code{lower}, and \code{upper} will also be returned, showing
#'   the information about the forecasts generated using all available
#'   observations.
#'
#' @family conformal prediction methods
#'
#' @seealso \code{\link{cvforecast}} to produce \code{object}, and
#'   \code{\link{update.cpforecast}} to extend the results with new
#'   observations.
#'
#' @references
#' Gibbs, I., and Candes, E. (2021). "Adaptive conformal inference under
#'   distribution shift", \emph{Advances in Neural Information Processing
#'   Systems}, \bold{34}, 1660--1672.
#'
#' @examples
#' # Simulate time series from an AR(2) model
#' library(forecast)
#' set.seed(1)
#' series <- arima.sim(n = 200, list(ar = c(0.8, -0.5)), sd = sqrt(1))
#'
#' # Cross-validation forecasting
#' far2 <- function(x, h, level) {
#'   Arima(x, order = c(2, 0, 0)) |>
#'     forecast(h = h, level)
#' }
#' fc <- cvforecast(series, forecastfun = far2, h = 3, level = 95, window = 50)
#'
#' # ACP with asymmetric nonconformity scores and rolling calibration sets
#' acpfc <- acp(fc, ncal = 50, rolling = TRUE)
#' print(acpfc)
#' summary(acpfc)
#'
#' @importFrom ggdist weighted_quantile
#' @export
acp <- function(
  object,
  alpha = 1 - 0.01 * object$level,
  gamma = 0.005,
  symmetric = FALSE,
  ncal = 10,
  rolling = FALSE,
  quantiletype = 1,
  update = FALSE,
  na.rm = TRUE,
  ...
) {
  # Check inputs
  if (any(alpha >= 1 | alpha <= 0)) {
    stop("alpha should be in (0, 1)")
  }
  if (gamma <= 0) {
    stop("the step size parameter gamma should be positive")
  }
  if (ncal < 10) {
    stop("length of calibration period should at least be 10")
  }
  if (!quantiletype %in% 1:9) {
    stop("quantiletype is invalid. It must be in 1:9.")
  }

  alpha <- sort(alpha, decreasing = TRUE)
  level <- 100 * (1 - alpha)
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

  namatrix <- ts(
    matrix(NA_real_, nrow = n, ncol = horizon),
    start = start(pf),
    frequency = frequency(pf)
  )
  colnames(namatrix) <- paste0("h=", seq(horizon))
  if (update) {
    lower <- object$LOWER
    upper <- object$UPPER
    if (symmetric) {
      alphat <- object$model$alpha_update$alpha
    } else {
      alphat_lower <- object$model$alpha_update$lower
      alphat_upper <- object$model$alpha_update$upper
    }
  } else {
    lower <- upper <-
      `names<-`(rep(list(namatrix), length(alpha)), paste0(level, "%"))
    if (symmetric) {
      alphat <- lower
    } else {
      alphat_lower <- alphat_upper <- lower
    }
  }

  out <- c(
    list(x = object$x, series = object$series),
    if ("xreg" %in% names(object)) list(xreg = object$xreg)
  )

  cp_times <- integer(horizon)
  for (h in seq(horizon)) {
    indx <- seq(ncal + h - 1, nrow(errors) - !object$forward, by = 1L)
    cp_times[h] <- length(indx)

    alphat_h <- alphat_lower_h <- alphat_upper_h <-
      errt_h <- errt_lower_h <- errt_upper_h <-
        q_lo_h <- q_up_h <-
          matrix(NA_real_, nrow = n, ncol = length(alpha))

    if (update) {
      for (i in seq(length(alpha))) {
        q_lo_h[, i] <- pf[, h] - lower[[i]][, h]
        q_up_h[, i] <- upper[[i]][, h] - pf[, h]
      }
      padded_errors <- rbind(
        errors,
        matrix(NA, nrow = n - nrow(errors), ncol = horizon)
      )

      if (symmetric) {
        alphat_h <- sapply(alphat, function(mat) mat[, h])
        errt_h <- abs(padded_errors[, h]) > q_lo_h
      } else {
        alphat_lower_h <- sapply(alphat_lower, function(mat) mat[, h])
        alphat_upper_h <- sapply(alphat_upper, function(mat) mat[, h])
        errt_lower_h <- (-padded_errors[, h]) > q_lo_h
        errt_upper_h <- padded_errors[, h] > q_up_h
      }
    }

    run_indx <- indx
    if (update) {
      rows <- indx + h
      missing_alpha <- if (symmetric) {
        rowSums(is.na(alphat_h[rows, , drop = FALSE])) > 0
      } else {
        rowSums(
          is.na(alphat_lower_h[rows, , drop = FALSE]) |
            is.na(alphat_upper_h[rows, , drop = FALSE])
        ) >
          0
      }
      first <- which(missing_alpha)[1L]
      if (is.na(first)) {
        run_indx <- integer()
      } else {
        first <- max(1L, first - 1L)
        run_indx <- indx[seq.int(first, length(indx))]
      }
    }

    for (t in run_indx) {
      errors_subset <- subset(
        errors[, h],
        start = ifelse(!rolling, 1, t - ncal + 1L),
        end = t
      )

      if (symmetric) {
        if (t == indx[1]) {
          alphat_h[t + h, ] <- alpha
        }

        # Compute sample quantiles
        q_lo_h[t + h, ] <- q_up_h[t + h, ] <- ggdist::weighted_quantile(
          x = abs(c(errors_subset, Inf)),
          probs = 1 - alphat_h[t + h, ],
          type = quantiletype,
          na.rm = na.rm,
          ...
        )

        if (t < tail(indx, 1)) {
          # Compute errt
          errt_h[t + 1, ] <- abs(errors[t + 1, h]) > q_lo_h[t + 1, ]
          outl <- which(alphat_h[t + 1, ] >= 1)
          outs <- which(alphat_h[t + 1, ] <= 0)
          errt_h[t + 1, outl] <- TRUE
          errt_h[t + 1, outs] <- FALSE

          if (anyNA(errt_h[t + 1, ])) {
            # Keep alpha unchanged
            alphat_h[t + h + 1, ] <- alphat_h[t + h, ]
          } else {
            # Update alpha
            alphat_h[t + h + 1, ] <- alphat_h[t + h, ] +
              gamma * (alpha - errt_h[t + 1, ])
          }
        }
      } else {
        if (t == indx[1]) {
          alphat_lower_h[t + h, ] <- alphat_upper_h[t + h, ] <- alpha / 2
        }

        # Compute sample quantiles
        q_lo_h[t + h, ] <- ggdist::weighted_quantile(
          x = -c(errors_subset, Inf),
          probs = 1 - alphat_lower_h[t + h, ],
          type = quantiletype,
          na.rm = na.rm,
          ...
        )
        q_up_h[t + h, ] <- ggdist::weighted_quantile(
          x = c(errors_subset, Inf),
          probs = 1 - alphat_upper_h[t + h, ],
          type = quantiletype,
          na.rm = na.rm,
          ...
        )

        if (t < tail(indx, 1)) {
          errt_lower_h[t + 1, ] <- (-errors[t + 1, h]) > q_lo_h[t + 1, ]
          errt_lower_h[t + 1, which(alphat_lower_h[t + 1, ] >= 1)] <- TRUE
          errt_lower_h[t + 1, which(alphat_lower_h[t + 1, ] <= 0)] <- FALSE

          errt_upper_h[t + 1, ] <- (errors[t + 1, h]) > q_up_h[t + 1, ]
          errt_upper_h[t + 1, which(alphat_upper_h[t + 1, ] >= 1)] <- TRUE
          errt_upper_h[t + 1, which(alphat_upper_h[t + 1, ] <= 0)] <- FALSE

          if (anyNA(errt_lower_h[t + 1, ]) || anyNA(errt_upper_h[t + 1, ])) {
            # Keep alpha unchanged
            alphat_lower_h[t + h + 1, ] <- alphat_lower_h[t + h, ]
            alphat_upper_h[t + h + 1, ] <- alphat_upper_h[t + h, ]
          } else {
            # Update alpha
            alphat_lower_h[t + h + 1, ] <- alphat_lower_h[t + h, ] +
              gamma * (alpha / 2 - errt_lower_h[t + 1, ])
            alphat_upper_h[t + h + 1, ] <- alphat_upper_h[t + h, ] +
              gamma * (alpha / 2 - errt_upper_h[t + 1, ])
          }
        }
      }
      for (i in seq(length(alpha))) {
        lower[[i]][t + h, h] <- pf[t + h, h] - q_lo_h[t + h, i]
        upper[[i]][t + h, h] <- pf[t + h, h] + q_up_h[t + h, i]
      }
    }
    for (i in seq(length(alpha))) {
      if (symmetric) {
        alphat[[i]][, h] <- alphat_h[, i]
      } else {
        alphat_lower[[i]][, h] <- alphat_lower_h[, i]
        alphat_upper[[i]][, h] <- alphat_upper_h[, i]
      }
    }
  }

  out$method <- paste("acp")
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
  out$model$args <- c(
    list(
      alpha = alpha,
      gamma = gamma,
      symmetric = symmetric,
      ncal = ncal,
      rolling = rolling,
      quantiletype = quantiletype,
      na.rm = na.rm
    ),
    list(...)
  )
  if (update) {
    out$model$cvforecast$call <- object$model$cvforecast$call
    out$model$cvforecast$args <- object$model$cvforecast$args
  } else {
    out$model$cvforecast$call <- object$call
    out$model$cvforecast$args <- object$args
  }
  if (symmetric) {
    out$model$alpha_update <- list(alpha = alphat)
  } else {
    out$model$alpha_update <- list(lower = alphat_lower, upper = alphat_upper)
  }

  return(structure(
    out,
    class = c("acp", "cpforecast", "cvforecast", "forecast")
  ))
}
