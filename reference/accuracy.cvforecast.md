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

- Winkler: Winkler Score

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
#>                ME       MAE      MSE      RMSE       MPE     MAPE      MASE
#> CV h=1 0.02360276 0.8586951 1.133966 0.8586951 110.01436 182.9608 0.7906662
#> CV h=2 0.05185524 1.1678584 2.006877 1.1678584  63.80171 188.6862 1.0730730
#> CV h=3 0.05616645 1.1728258 2.021815 1.1728258  46.44065 181.8183 1.0767439
#>            RMSSE
#> CV h=1 0.6306113
#> CV h=2 0.8562614
#> CV h=3 0.8593417
accuracy(fc, measures = interval_measures, level = 95, byhorizon = TRUE)
#>        Winkler_95  MSIS_95
#> CV h=1   5.137519 4.708273
#> CV h=2   5.992869 5.493790
#> CV h=3   6.157379 5.646547

# Out-of-sample forecast accuracy on test set
accuracy(fc, x = c(1, 0.5, 0), measures = interval_measures,
         level = 95, byhorizon = TRUE)
#>          Winkler_95  MSIS_95
#> CV h=1     5.137519 4.708273
#> CV h=2     5.992869 5.493790
#> CV h=3     6.157379 5.646547
#> Test h=1   4.723912 4.386664
#> Test h=2   5.809702 5.394939
#> Test h=3   5.841090 5.424085
```
