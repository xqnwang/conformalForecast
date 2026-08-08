# Conformal prediction

This function allows you to specify the method used to perform conformal
prediction.

## Usage

``` r
conformal(object, ...)

# S3 method for class 'cvforecast'
conformal(object, method = c("scp", "acp", "pid", "acmcp"), ...)
```

## Arguments

- object:

  An object of class `"cvforecast"`. It must have an argument `x` for
  the original univariate time series, an argument `MEAN` for the point
  forecasts and `ERROR` for the forecast errors on the validation set.
  See the results of a call to
  [`cvforecast`](https://xqnwang.github.io/conformalForecast/reference/cvforecast.md).

- ...:

  Additional arguments to be passed to the selected conformal method.

- method:

  A character string specifying the conformal method to be applied.
  Possible options include `"scp"`
  ([scp](https://xqnwang.github.io/conformalForecast/reference/scp.md)),
  `"acp"`
  ([acp](https://xqnwang.github.io/conformalForecast/reference/acp.md)),
  `"pid"`
  ([pid](https://xqnwang.github.io/conformalForecast/reference/pid.md)),
  and `"acmcp"`
  ([acmcp](https://xqnwang.github.io/conformalForecast/reference/acmcp.md)).
  Defaults to `"scp"`.

## Value

An object whose class depends on the method invoked.

## See also

[`cvforecast`](https://xqnwang.github.io/conformalForecast/reference/cvforecast.md)
to produce `object`, and
[`update.cpforecast`](https://xqnwang.github.io/conformalForecast/reference/update.cpforecast.md)
to extend the results with new observations.

Other conformal prediction methods:
[`acmcp()`](https://xqnwang.github.io/conformalForecast/reference/acmcp.md),
[`acp()`](https://xqnwang.github.io/conformalForecast/reference/acp.md),
[`pid()`](https://xqnwang.github.io/conformalForecast/reference/pid.md),
[`scp()`](https://xqnwang.github.io/conformalForecast/reference/scp.md)

## Examples

``` r
# Simulate time series from an AR(2) model
library(forecast)
set.seed(1)
series <- arima.sim(n = 200, list(ar = c(0.8, -0.5)), sd = sqrt(1))

# Cross-validation forecasting
far2 <- function(x, h, level) {
  Arima(x, order = c(2, 0, 0)) |>
    forecast(h = h, level)
}
fc <- cvforecast(series, forecastfun = far2, h = 3, level = 95, window = 50)

# Classical conformal prediction with equal weights
scpfc <- conformal(fc, method = "scp", ncal = 50, rolling = TRUE)
summary(scpfc)
#> SCP 
#> 
#> Call:
#>  scp(object = object, ncal = 50, rolling = TRUE) 
#> 
#>  cp_times (the forward step included): 101 (h=1), 100 (h=2), 99 (h=3)
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201      0.6271538 -1.023525 3.234253
#> 202      0.8607034 -1.291456 4.064625
#> 203      0.4935805 -1.675089 3.691253
#> 
#> Cross-validation error measures:
#>       ME   MAE   MSE RMSE    MPE    MAPE  MASE RMSSE Winkler_95 MSIS_95
#> CV 0.007 0.946 1.415 1.06 -3.933 269.763 0.992 0.882      6.123   6.568

# ACP with asymmetric nonconformity scores and rolling calibration sets
acpfc <- conformal(fc, method = "acp", gamma = 0.005,
                   ncal = 50, rolling = TRUE)
summary(acpfc)
#> ACP 
#> 
#> Call:
#>  acp(object = object, gamma = 0.005, ncal = 50, rolling = TRUE) 
#> 
#>  cp_times (the forward step included): 101 (h=1), 100 (h=2), 99 (h=3)
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201      0.6271538 -1.706076 3.234253
#> 202      0.8607034 -1.825564      Inf
#> 203      0.4935805 -2.027163 3.691253
#> 
#> Cross-validation error measures:
#>       ME   MAE   MSE RMSE    MPE    MAPE  MASE RMSSE Winkler_95 MSIS_95
#> CV 0.007 0.946 1.415 1.06 -3.933 269.763 0.992 0.882        Inf     Inf
```
