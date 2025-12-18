# Accuracy measures for a cross-validation model and a conformal prediction model

Return range of summary measures of the out-of-sample forecast accuracy.
If `x` is given, the function also measures test set forecast accuracy.
If `x` is not given, the function only produces accuracy measures on
validation set.

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

  An optional numerical vector containing actual values of the same
  length as `mean` in `object`.

- CV:

  If `TRUE`, the cross-validation forecast accuracy will be returned.

- period:

  The seasonal period of the data.

- measures:

  A list of accuracy measure functions to compute (such as
  [point_measures](https://xqnwang.github.io/conformalForecast/reference/point_measures.md)
  or
  [interval_measures](https://xqnwang.github.io/conformalForecast/reference/interval_measures.md)).

- byhorizon:

  If `TRUE`, accuracy measures will be calculated for each individual
  forecast horizon `h` separately.

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

- winkler_score: Winkler Score

- MSIS: Mean Scaled Interval Score

## See also

[`point_measures`](https://xqnwang.github.io/conformalForecast/reference/point_measures.md),
[`interval_measures`](https://xqnwang.github.io/conformalForecast/reference/interval_measures.md)

## Examples

``` r
# Simulate time series from an AR(2) model
library(forecast)
series <- arima.sim(n = 200, list(ar = c(0.8, -0.5)), sd = sqrt(1))

# Cross-validation forecasting with a rolling window
far2 <- function(x, h, level) {
  Arima(x, order = c(2, 0, 0)) |>
    forecast(h = h, level)
}
fc <- cvforecast(series, forecastfun = far2, h = 3, level = 95,
                 forward = TRUE, initial = 1, window = 50)

# Out-of-sample forecast accuracy on validation set
accuracy(fc, measures = point_measures, byhorizon = TRUE)
#>                ME       MAE      MSE      RMSE       MPE     MAPE     MASE
#> CV h=1 0.03097937 0.8513179 1.119499 0.8513179 111.45557 184.4020 0.778602
#> CV h=2 0.05953212 1.1602238 1.998144 1.1602238  63.31647 188.2276 1.058245
#> CV h=3 0.04710356 1.1675812 2.014782 1.1675812  47.17800 182.5556 1.063932
#>            RMSSE
#> CV h=1 0.6228594
#> CV h=2 0.8471084
#> CV h=3 0.8518445
accuracy(fc, measures = interval_measures, level = 95, byhorizon = TRUE)
#>        Winkler_95  MSIS_95
#> CV h=1   5.141322 4.677845
#> CV h=2   5.997589 5.458960
#> CV h=3   6.161738 5.610096

# Out-of-sample forecast accuracy on test set
accuracy(fc, x = c(1, 0.5, 0), measures = interval_measures,
         level = 95, byhorizon = TRUE)
#>          Winkler_95  MSIS_95
#> CV h=1     5.141322 4.677845
#> CV h=2     5.997589 5.458960
#> CV h=3     6.161738 5.610096
#> Test h=1   4.678008 4.341244
#> Test h=2   5.807041 5.389000
#> Test h=3   5.850943 5.429741
```
