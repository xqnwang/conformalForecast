# Package index

## Forecasting workflow

Generate multistep-ahead forecasts by time series cross-validation,
conformalize them, and update the result as new observations arrive.

- [`cvforecast()`](https://xqnwang.github.io/conformalForecast/reference/cvforecast.md)
  : Time series cross-validation forecasting
- [`conformal()`](https://xqnwang.github.io/conformalForecast/reference/conformal.md)
  : Conformal prediction
- [`update(`*`<cpforecast>`*`)`](https://xqnwang.github.io/conformalForecast/reference/update.cpforecast.md)
  : Update and reperform cross-validation forecasting and conformal
  prediction

## Conformal prediction methods

The calibration algorithms dispatched by
[`conformal()`](https://xqnwang.github.io/conformalForecast/reference/conformal.md).

- [`scp()`](https://xqnwang.github.io/conformalForecast/reference/scp.md)
  : Classical split conformal prediction method
- [`acp()`](https://xqnwang.github.io/conformalForecast/reference/acp.md)
  : Adaptive conformal prediction method
- [`pid()`](https://xqnwang.github.io/conformalForecast/reference/pid.md)
  : Conformal PID control method
- [`acmcp()`](https://xqnwang.github.io/conformalForecast/reference/acmcp.md)
  : Autocorrelated multistep-ahead conformal prediction method

## Forecast evaluation

Assess point and interval forecast performance.

- [`accuracy(`*`<cvforecast>`*`)`](https://xqnwang.github.io/conformalForecast/reference/accuracy.cvforecast.md)
  [`accuracy(`*`<cpforecast>`*`)`](https://xqnwang.github.io/conformalForecast/reference/accuracy.cvforecast.md)
  : Accuracy measures for a cross-validation model and a conformal
  prediction model
- [`coverage()`](https://xqnwang.github.io/conformalForecast/reference/coverage.md)
  : Calculate interval forecast coverage
- [`width()`](https://xqnwang.github.io/conformalForecast/reference/width.md)
  : Calculate interval forecast width
- [`ME()`](https://xqnwang.github.io/conformalForecast/reference/point_measures.md)
  [`MAE()`](https://xqnwang.github.io/conformalForecast/reference/point_measures.md)
  [`MSE()`](https://xqnwang.github.io/conformalForecast/reference/point_measures.md)
  [`RMSE()`](https://xqnwang.github.io/conformalForecast/reference/point_measures.md)
  [`MPE()`](https://xqnwang.github.io/conformalForecast/reference/point_measures.md)
  [`MAPE()`](https://xqnwang.github.io/conformalForecast/reference/point_measures.md)
  [`MASE()`](https://xqnwang.github.io/conformalForecast/reference/point_measures.md)
  [`RMSSE()`](https://xqnwang.github.io/conformalForecast/reference/point_measures.md)
  [`point_measures`](https://xqnwang.github.io/conformalForecast/reference/point_measures.md)
  : Point estimate accuracy measures
- [`MSIS()`](https://xqnwang.github.io/conformalForecast/reference/interval_measures.md)
  [`winkler_score()`](https://xqnwang.github.io/conformalForecast/reference/interval_measures.md)
  [`interval_measures`](https://xqnwang.github.io/conformalForecast/reference/interval_measures.md)
  : Interval estimate accuracy measures

## Utilities

- [`lagmatrix()`](https://xqnwang.github.io/conformalForecast/reference/lagmatrix.md)
  : Create lags or leads of a matrix

## Package overview

- [`conformalForecast`](https://xqnwang.github.io/conformalForecast/reference/conformalForecast-package.md)
  [`conformalForecast-package`](https://xqnwang.github.io/conformalForecast/reference/conformalForecast-package.md)
  : conformalForecast: Conformal Prediction Methods for Multistep-Ahead
  Time Series Forecasting
