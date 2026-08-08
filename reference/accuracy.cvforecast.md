# Accuracy measures for a cross-validation model and a conformal prediction model

Return a range of summary measures of the out-of-sample forecast
accuracy. If `x` is given, the function also measures the test set
forecast accuracy. If `x` is not given, the function only produces
accuracy measures on the validation set.

## Usage

``` r
# S3 method for class 'cvforecast'
accuracy(
  object,
  x,
  CV = TRUE,
  period = NULL,
  measures = interval_measures,
  byhorizon = FALSE,
  ...
)

# S3 method for class 'cpforecast'
accuracy(object, ...)
```

## Arguments

- object:

  An object of class `"cvforecast"` or `"cpforecast"`.

- x:

  An optional numerical vector containing the actual values of the same
  length as `mean` in `object`.

- CV:

  If `TRUE`, the cross-validation forecast accuracy will be returned.
  Defaults to `TRUE`.

- period:

  The seasonal period of the data. Defaults to `NULL`, in which case it
  is taken from the frequency of the series.

- measures:

  A list of accuracy measure functions to compute (such as
  [point_measures](https://xqnwang.github.io/conformalForecast/reference/point_measures.md)
  or
  [interval_measures](https://xqnwang.github.io/conformalForecast/reference/interval_measures.md)).
  Defaults to `interval_measures`.

- byhorizon:

  If `TRUE`, accuracy measures will be calculated for each individual
  forecast horizon `h` separately. Defaults to `FALSE`.

- ...:

  Additional arguments depending on the specific measure.

## Value

A matrix giving mean out-of-sample forecast accuracy measures.

## Details

The measures calculated are:

- ME: Mean Error

- MAE: Mean Absolute Error

- MSE: Mean Squared Error

- RMSE: Root Mean Squared Error

- MPE: Mean Percentage Error

- MAPE: Mean Absolute Percentage Error

- MASE: Mean Absolute Scaled Error

- RMSSE: Root Mean Squared Scaled Error

- Winkler: Winkler Score

- MSIS: Mean Scaled Interval Score

## See also

[`point_measures`](https://xqnwang.github.io/conformalForecast/reference/point_measures.md)
and
[`interval_measures`](https://xqnwang.github.io/conformalForecast/reference/interval_measures.md)
for the measures that can be supplied through `measures`.

Other evaluation functions:
[`coverage()`](https://xqnwang.github.io/conformalForecast/reference/coverage.md),
[`width()`](https://xqnwang.github.io/conformalForecast/reference/width.md)

## Examples

``` r
# Simulate time series from an AR(2) model
library(forecast)
set.seed(1)
series <- arima.sim(n = 200, list(ar = c(0.8, -0.5)), sd = sqrt(1))

# Cross-validation forecasting with a rolling window
far2 <- function(x, h, level) {
  Arima(x, order = c(2, 0, 0)) |>
    forecast(h = h, level)
}
fc <- cvforecast(series, forecastfun = far2, h = 3, level = 95, window = 50)

# Out-of-sample forecast accuracy on validation set
accuracy(fc, measures = point_measures, byhorizon = TRUE)
#>                 ME       MAE      MSE      RMSE       MPE     MAPE     MASE
#> CV h=1  0.02185523 0.8052604 1.017906 0.8052604  35.66722 221.8293 0.844289
#> CV h=2  0.01698001 1.0245970 1.632225 1.0245970 -32.41324 319.0622 1.075412
#> CV h=3 -0.01246597 1.0183022 1.622138 1.0183022 -17.90026 271.4276 1.067876
#>            RMSSE
#> CV h=1 0.6697150
#> CV h=2 0.8532125
#> CV h=3 0.8477134
accuracy(fc, measures = interval_measures, level = 95, byhorizon = TRUE)
#>        Winkler_95  MSIS_95
#> CV h=1   5.012254 5.296003
#> CV h=2   6.407682 6.768502
#> CV h=3   6.245013 6.565268

# Accuracy of conformal prediction intervals
scpfc <- conformal(fc, method = "scp", ncal = 50, rolling = TRUE)
accuracy(scpfc, measures = interval_measures, level = 95,
         byhorizon = TRUE)
#>        Winkler_95  MSIS_95
#> CV h=1   5.394495 5.789155
#> CV h=2   6.629521 7.129421
#> CV h=3   6.502841 6.963483

# Out-of-sample forecast accuracy on test set
accuracy(fc, x = c(1, 0.5, 0), measures = interval_measures,
         CV = FALSE, level = 95, byhorizon = TRUE)
#>          Winkler_95  MSIS_95
#> Test h=1   4.079088 4.073511
#> Test h=2   5.045214 5.038316
#> Test h=3   5.050806 5.043901
```
