# conformalForecast (development version)

* Settings previously stored as `model$alpha`, `model$symmetric`,
  `model$integrate`, `model$scorecast`, `model$lr`, `model$Csat` and `model$KI`
  are now collected in `model$args`.
* `forecast` moved from `Imports` to `Depends`, so it is attached together with
  `conformalForecast`.
* `acmcp()` now accepts `ma_method = "CSS"` for faster MA scorecaster fitting at
  longer forecast horizons, while retaining `"CSS-ML"` as the default, and
  reports successful scorecasts in `scorecast_times`.
* New unit-test suite covering all four conformal methods, `cvforecast()`, 
  `update()`, and the evaluation functions.

## Bug fixes

* `acmcp()` no longer fits MA models before its recursive scorecaster inputs are
  available.
* `update()` now works on `pid()`, `acmcp()` and `acp(symmetric = TRUE)`
  objects, and resumes from the last computed step instead of recomputing the
  whole history.
* `update()` now replays stored arguments, preserves forecasts at all
  frequencies, validates new data and external regressors, and handles failed
  final model fits. External regressors are retained by `pid()` and `acmcp()`,
  enabling external-regressor updates for all four conformal methods.
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
* Expanded the vignette with the `conformal()` interface, incremental updates,
  external regressors, clearer method descriptions, and updated ggplot2 calls.
* New `vignette("update")`, on extending a fitted conformal pipeline with newly
  arrived observations instead of refitting it from scratch.
* The introductory vignette now documents `conformal()`, and describes each
  method, its key formula and its tuning arguments in more detail.
* The pkgdown site indexes both vignettes under "Articles".
* Help examples now cover raw inputs to `coverage()` and `width()` and
  `cpforecast` accuracy, and use fixed random seeds for reproducibility.

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
