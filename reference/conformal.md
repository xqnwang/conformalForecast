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
  original univariate time series, an argument `MEAN` for point
  forecasts and `ERROR` for forecast errors on validation set. See the
  results of a call to
  [`cvforecast`](https://xqnwang.github.io/conformalForecast/reference/cvforecast.md).

- ...:

  Additional arguments to be passed to the selected conformal method.

- method:

  A character string specifying the conformal method to be applied.
  Possible options include `"scp"`
  ([scp](https://xqnwang.github.io/conformalForecast/reference/scp.md)),
  `"acp"`([acp](https://xqnwang.github.io/conformalForecast/reference/acp.md)),
  `"pid"`([pid](https://xqnwang.github.io/conformalForecast/reference/pid.md)),
  and
  `"acmcp"`([acmcp](https://xqnwang.github.io/conformalForecast/reference/acmcp.md)).

## Value

An object whose class depends on the method invoked.

## See also

[`scp`](https://xqnwang.github.io/conformalForecast/reference/scp.md),
[`acp`](https://xqnwang.github.io/conformalForecast/reference/acp.md),
[`pid`](https://xqnwang.github.io/conformalForecast/reference/pid.md),
and
[`acmcp`](https://xqnwang.github.io/conformalForecast/reference/acmcp.md)

## Examples

``` r
# Simulate time series from an AR(2) model
library(forecast)
series <- arima.sim(n = 200, list(ar = c(0.8, -0.5)), sd = sqrt(1))

# Cross-validation forecasting
far2 <- function(x, h, level) {
  Arima(x, order = c(2, 0, 0)) |>
    forecast(h = h, level)
}
fc <- cvforecast(series, forecastfun = far2, h = 3, level = 95,
                 forward = TRUE, initial = 1, window = 50)

# Classical conformal prediction with equal weights
scpfc <- conformal(fc, method = "scp", symmetric = FALSE, ncal = 50, rolling = TRUE)
summary(scpfc)
#> SCP 
#> 
#> Call:
#>  scp(object = object, symmetric = FALSE, ncal = 50, rolling = TRUE) 
#> 
#>  cp_times = 99 (the forward step included) 
#> 
#> Forecasts of the forward step:
#> Cross-validation
#> 
#> Call:
#>  scp(object = object, symmetric = FALSE, ncal = 50, rolling = TRUE) 
#> 
#>  fit_times =  (the forward step included) 
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201      0.3481845 -1.694398 2.606816
#> 202     -0.4556896 -2.590600 1.370540
#> 203     -0.7068473 -3.113914 1.136993
#> 
#> Cross-validation error measures:
#>        ME   MAE   MSE RMSE     MPE    MAPE  MASE RMSSE Winkler_95 MSIS_95
#> CV -0.001 0.994 1.501  1.1 177.356 266.543 0.902 0.784      7.043   6.409

# ACP with asymmetric nonconformity scores and rolling calibration sets
acpfc <- conformal(fc, method = "acp", symmetric = FALSE, gamma = 0.005,
                   ncal = 50, rolling = TRUE)
summary(acpfc)
#> ACP 
#> 
#> Call:
#>  acp(object = object, gamma = 0.005, symmetric = FALSE, ncal = 50,  
#>      rolling = TRUE) 
#> 
#>  cp_times = 99 (the forward step included) 
#> 
#> Forecasts of the forward step:
#> Cross-validation
#> 
#> Call:
#>  acp(object = object, gamma = 0.005, symmetric = FALSE, ncal = 50,  
#>      rolling = TRUE) 
#> 
#>  fit_times =  (the forward step included) 
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201      0.3481845 -1.694398 2.606816
#> 202     -0.4556896 -4.358829 1.370540
#> 203     -0.7068473 -4.577801      Inf
#> 
#> Cross-validation error measures:
#>        ME   MAE   MSE RMSE     MPE    MAPE  MASE RMSSE Winkler_95 MSIS_95
#> CV -0.001 0.994 1.501  1.1 177.356 266.543 0.902 0.784        Inf     Inf
```
