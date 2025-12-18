# Autocorrelated multistep-ahead conformal prediction method

Compute prediction intervals and other information by applying the
Autocorrelated Multistep-ahead Conformal Prediction (AcMCP) method. The
method can only deal with asymmetric nonconformity scores, i.e.,
forecast errors.

## Usage

``` r
acmcp(
  object,
  alpha = 1 - 0.01 * object$level,
  ncal = 10,
  rolling = FALSE,
  integrate = TRUE,
  scorecast = TRUE,
  lr = 0.1,
  Tg = NULL,
  delta = NULL,
  Csat = 2/pi * (ceiling(log(Tg) * delta) - 1/log(Tg)),
  KI = max(abs(object$errors), na.rm = TRUE),
  ...
)
```

## Arguments

- object:

  An object of class `"cvforecast"`. It must have an argument `x` for
  original univariate time series, an argument `MEAN` for point
  forecasts and `ERROR` for forecast errors on validation set. See the
  results of a call to
  [`cvforecast`](https://xqnwang.github.io/conformalForecast/reference/cvforecast.md).

- alpha:

  A numeric vector of significance levels to achieve a desired coverage
  level \\1-\alpha\\.

- ncal:

  Length of the burn-in period for training the scorecaster. If
  `rolling = TRUE`, it is also used as the length of the trailing
  windows for learning rate calculation and the windows for the
  calibration set. If `rolling = FALSE`, it is used as the initial
  period of calibration sets and trailing windows for learning rate
  calculation.

- rolling:

  If `TRUE`, a rolling window strategy will be adopted to form the
  trailing window for learning rate calculation and the calibration set
  for scorecaster if applicable. Otherwise, expanding window strategy
  will be used.

- integrate:

  If `TRUE`, error integration will be included in the update process.

- scorecast:

  If `TRUE`, scorecasting will be included in the update process.

- lr:

  Initial learning rate used for quantile tracking.

- Tg:

  The time that is set to achieve the target absolute coverage guarantee
  before this.

- delta:

  The target absolute coverage guarantee is set to \\1-\alpha-\delta\\.

- Csat:

  A positive constant ensuring that by time `Tg`, an absolute guarantee
  is of at least \\1-\alpha-\delta\\ coverage.

- KI:

  A positive constant to place the integrator on the same scale as the
  scores.

- ...:

  Other arguments are passed to the function.

## Value

A list of class `c("acmcp", "cpforecast", "forecast")` with the
following components:

- x:

  The original time series.

- series:

  The name of the series `x`.

- method:

  A character string "acmcp".

- cp_times:

  The number of times the conformal prediction is performed in
  cross-validation.

- MEAN:

  Point forecasts as a multivariate time series, where the \\h\\th
  column holds the point forecasts for forecast horizon \\h\\. The time
  index corresponds to the period for which the forecast is produced.

- ERROR:

  Forecast errors given by \\e\_{t+h\|t} = y\_{t+h}-\hat{y}\_{t+h\|t}\\.

- LOWER:

  A list containing lower bounds for prediction intervals for each
  `level`. Each element within the list will be a multivariate time
  series with the same dimensional characteristics as `MEAN`.

- UPPER:

  A list containing upper bounds for prediction intervals for each
  `level`. Each element within the list will be a multivariate time
  series with the same dimensional characteristics as `MEAN`.

- level:

  The confidence values associated with the prediction intervals.

- call:

  The matched call.

- model:

  A list containing information abouth the conformal prediction model.

If `mean` is included in the `object`, the components `mean`, `lower`,
and `upper` will also be returned, showing the information about the
test set forecasts generated using all available observations.

## Details

Similar to the PID method, the AcMCP method also integrates three
modules (P, I, and D) to form the final iteration. However, instead of
performing conformal prediction for each individual forecast horizon `h`
separately, AcMCP employs a combination of an MA\\(h-1)\\ model and a
linear regression model of \\e\_{t+h\|t}\\ on
\\e\_{t+h-1\|t},\dots,e\_{t+1\|t}\\ as the scorecaster. This allows the
AcMCP method to capture the relationship between the \\h\\-step ahead
forecast error and past errors.

## References

Wang, X., and Hyndman, R. J. (2024). "Online conformal inference for
multi-step time series forecasting", arXiv preprint arXiv:2410.13115.

## See also

[`pid`](https://xqnwang.github.io/conformalForecast/reference/pid.md)

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

# AcMCP setup
Tg <- 200; delta <- 0.01
Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
KI <- 2
lr <- 0.1

# AcMCP with integrator and scorecaster
acmcpfc <- acmcp(fc, ncal = 50, rolling = TRUE,
             integrate = TRUE, scorecast = TRUE,
             lr = lr, KI = KI, Csat = Csat)
print(acmcpfc)
#> ACMCP 
#> 
#> Call:
#>  acmcp(object = fc, ncal = 50, rolling = TRUE, integrate = TRUE,  
#>      scorecast = TRUE, lr = lr, Csat = Csat, KI = KI) 
#> 
#>  cp_times = 148 (the forward step included) 
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201      0.5094989 -1.732309 2.155511
#> 202      0.4044399 -3.233107 2.002388
#> 203     -0.1040782 -3.529078 1.363847
summary(acmcpfc)
#> ACMCP 
#> 
#> Call:
#>  acmcp(object = fc, ncal = 50, rolling = TRUE, integrate = TRUE,  
#>      scorecast = TRUE, lr = lr, Csat = Csat, KI = KI) 
#> 
#>  cp_times = 148 (the forward step included) 
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201      0.5094989 -1.732309 2.155511
#> 202      0.4044399 -3.233107 2.002388
#> 203     -0.1040782 -3.529078 1.363847
#> 
#> Cross-validation error measures:
#>        ME   MAE   MSE  RMSE   MPE    MAPE  MASE RMSSE Winkler_95 MSIS_95
#> CV -0.109 1.003 1.566 1.143 -4.17 689.961 0.905 0.826      6.258   5.532
```
