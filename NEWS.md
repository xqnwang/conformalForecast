# conformalForecast (development version)

* Settings previously stored as `model$alpha`, `model$symmetric`,
  `model$integrate`, `model$scorecast`, `model$lr`, `model$Csat` and `model$KI`
  are now collected in `model$args`.
* `forecast` moved from `Imports` to `Depends`, so it is attached together with
  `conformalForecast`.

## Bug fixes

* `update()` now works on `pid()`, `acmcp()` and `acp(symmetric = TRUE)`
  objects, and is much faster: it resumes from the last computed step instead of
  recomputing the whole history.
* `update()` no longer depends on the calling environment, and no longer
  silently drops the new forecasts at some frequencies.
* `pid()` and `acmcp()` can be called with their documented defaults.
* `print()` and `summary()` no longer repeat the cross-validation header.
* `cvforecast()` no longer warns about arguments legitimately passed to
  `forecastfun`, and survives a failing final model fit.
* `cp_times` is reported per forecast horizon.
* `coverage()` and `width()` accept `x`, `LOWER` and `UPPER` through `...`, and
  require a single `level`.
* `lagmatrix()` handles lags larger than the number of rows.
* Clearer error from `scp()` and `acp()` when `ncal` is too large.

## Documentation

* The pkgdown site renders mathematics again (`math-rendering: mathjax`).
* Corrected the documented return classes, the `Winkler` measure name, the
  `lagmatrix()` error message, and several typos.
* Vignette: replaced deprecated ggplot2 calls.
* Fixed several typos.

# conformalForecast 0.1.1

* Changed accuracy.default to use S3 methods

# conformalForecast 0.1.0

* First release.

## Initial features

* Provides implementations of the main conformal prediction methods:
  - `scp` (split conformal prediction)
  - `acp` (adaptive conformal prediction)
  - `pid` (conformal PID control method)
  - `acmcp` (autocorrelated multi-step conformal prediction).
* Includes utility functions to compute coverage rates and interval widths for evaluating predictive performance.
