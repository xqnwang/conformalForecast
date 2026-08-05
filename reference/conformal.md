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
#>  cp_times (the forward step included): 101 (h=1), 100 (h=2), 99 (h=3)
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201    -0.80664167 -2.849224 1.451989
#> 202    -0.20322450 -2.338135 1.623005
#> 203     0.05669504 -2.350372 1.900536
#> 
#> Cross-validation error measures:
#>        ME   MAE   MSE  RMSE     MPE    MAPE MASE RMSSE Winkler_95 MSIS_95
#> CV -0.042 0.997 1.523 1.105 175.173 265.181 0.89  0.78      7.076   6.359

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
#>  cp_times (the forward step included): 101 (h=1), 100 (h=2), 99 (h=3)
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201    -0.80664167 -2.849224 1.451989
#> 202    -0.20322450 -4.106364 1.623005
#> 203     0.05669504 -3.814259      Inf
#> 
#> Cross-validation error measures:
#>        ME   MAE   MSE  RMSE     MPE    MAPE MASE RMSSE Winkler_95 MSIS_95
#> CV -0.042 0.997 1.523 1.105 175.173 265.181 0.89  0.78        Inf     Inf
```
