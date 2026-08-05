# Changelog

## conformalForecast (development version)

- Settings previously stored as `model$alpha`, `model$symmetric`,
  `model$integrate`, `model$scorecast`, `model$lr`, `model$Csat` and
  `model$KI` are now collected in `model$args`.
- `forecast` moved from `Imports` to `Depends`, so it is attached
  together with `conformalForecast`.

### Bug fixes

- [`update()`](https://rdrr.io/r/stats/update.html) now works on
  [`pid()`](https://xqnwang.github.io/conformalForecast/reference/pid.md),
  [`acmcp()`](https://xqnwang.github.io/conformalForecast/reference/acmcp.md)
  and `acp(symmetric = TRUE)` objects, and resumes from the last
  computed step instead of recomputing the whole history.
- [`update()`](https://rdrr.io/r/stats/update.html) now replays stored
  arguments, preserves forecasts at all frequencies, validates new data
  and external regressors, and handles failed final model fits. External
  regressors are retained by
  [`pid()`](https://xqnwang.github.io/conformalForecast/reference/pid.md)
  and
  [`acmcp()`](https://xqnwang.github.io/conformalForecast/reference/acmcp.md),
  enabling external-regressor updates for all four conformal methods.
- [`pid()`](https://xqnwang.github.io/conformalForecast/reference/pid.md)
  and
  [`acmcp()`](https://xqnwang.github.io/conformalForecast/reference/acmcp.md)
  can be called with their documented defaults.
- [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) no longer repeat
  the cross-validation header.
- [`cvforecast()`](https://xqnwang.github.io/conformalForecast/reference/cvforecast.md)
  no longer warns about arguments legitimately passed to `forecastfun`,
  and survives a failing final model fit.
- `cp_times` is reported per forecast horizon.
- [`coverage()`](https://xqnwang.github.io/conformalForecast/reference/coverage.md)
  and
  [`width()`](https://xqnwang.github.io/conformalForecast/reference/width.md)
  accept `x`, `LOWER` and `UPPER` through `...`, and require a single
  `level`.
- [`lagmatrix()`](https://xqnwang.github.io/conformalForecast/reference/lagmatrix.md)
  handles lags larger than the number of rows.
- Clearer error from
  [`scp()`](https://xqnwang.github.io/conformalForecast/reference/scp.md)
  and
  [`acp()`](https://xqnwang.github.io/conformalForecast/reference/acp.md)
  when `ncal` is too large.

### Documentation

- The pkgdown site renders mathematics again
  (`math-rendering: mathjax`).
- Corrected the documented return classes, the `Winkler` measure name,
  the
  [`lagmatrix()`](https://xqnwang.github.io/conformalForecast/reference/lagmatrix.md)
  error message, and several typos.
- Expanded the vignette with the
  [`conformal()`](https://xqnwang.github.io/conformalForecast/reference/conformal.md)
  interface, incremental updates, external regressors, clearer method
  descriptions, and updated ggplot2 calls.
- Help examples now cover raw inputs to
  [`coverage()`](https://xqnwang.github.io/conformalForecast/reference/coverage.md)
  and
  [`width()`](https://xqnwang.github.io/conformalForecast/reference/width.md)
  and `cpforecast` accuracy, and use fixed random seeds for
  reproducibility.

## conformalForecast 0.1.1

CRAN release: 2026-01-15

- Changed accuracy.default to use S3 methods

## conformalForecast 0.1.0

CRAN release: 2025-10-06

- First release.

### Initial features

- Provides implementations of the main conformal prediction methods:
  - `scp` (split conformal prediction)
  - `acp` (adaptive conformal prediction)
  - `pid` (conformal PID control method)
  - `acmcp` (autocorrelated multi-step conformal prediction).
- Includes utility functions to compute coverage rates and interval
  widths for evaluating predictive performance.
