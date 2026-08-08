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
  Tg = NROW(object$ERROR),
  delta = 0.01,
  Csat = NULL,
  KI = max(abs(object$ERROR), na.rm = TRUE),
  update = FALSE,
  ma_method = c("CSS-ML", "CSS"),
  ...
)
```

## Arguments

- object:

  An object of class `"cvforecast"`. It must have an argument `x` for
  the original univariate time series, an argument `MEAN` for the point
  forecasts and `ERROR` for the forecast errors on the validation set.
  See the results of a call to
  [`cvforecast`](https://xqnwang.github.io/conformalForecast/reference/cvforecast.md).

- alpha:

  A numeric vector of significance levels to achieve a desired coverage
  level \\1-\alpha\\. Defaults to `1 - 0.01 * object$level`, the levels
  used in `object`.

- ncal:

  Length of the burn-in period for training the scorecaster. If
  `rolling = TRUE`, it is also used as the length of the trailing
  windows for the learning rate calculation and of the windows for the
  calibration set. If `rolling = FALSE`, it is used as the initial
  period of the calibration sets and of the trailing windows for the
  learning rate calculation. Defaults to `10`.

- rolling:

  If `TRUE`, a rolling window strategy will be adopted to form the
  trailing window for the learning rate calculation and the calibration
  set for the scorecaster if applicable. Otherwise, an expanding window
  strategy will be used. Defaults to `FALSE`.

- integrate:

  If `TRUE`, error integration will be included in the update process.
  Defaults to `TRUE`.

- scorecast:

  If `TRUE`, scorecasting will be included in the update process.
  Defaults to `TRUE`.

- lr:

  A positive initial learning rate used for quantile tracking. Defaults
  to `0.1`.

- Tg:

  The time that is set to achieve the target absolute coverage guarantee
  before this. It must be greater than 1 when `Csat` is not supplied.
  Defaults to `NROW(object$ERROR)`, the number of cross-validation
  periods in `object`.

- delta:

  A number in \\(0, 1)\\. The target absolute coverage guarantee is set
  to \\1-\alpha-\delta\\. Defaults to `0.01`.

- Csat:

  A positive constant ensuring that by time `Tg`, an absolute guarantee
  is of at least \\1-\alpha-\delta\\ coverage. Defaults to `NULL`, in
  which case it is derived from `Tg` and `delta` as
  `2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))`.

- KI:

  A non-negative constant to place the integrator on the same scale as
  the scores. Defaults to `max(abs(object$ERROR), na.rm = TRUE)`, the
  largest absolute forecast error in `object`.

- update:

  If `TRUE`, `object` already holds the results of a previous call and
  only the newly added time steps are computed; the prediction intervals
  produced earlier are carried over unchanged. Set by
  [`update.cpforecast`](https://xqnwang.github.io/conformalForecast/reference/update.cpforecast.md)
  and not normally set by hand. Defaults to `FALSE`.

- ma_method:

  Estimation method for the MA\\(h-1)\\ scorecaster. `"CSS-ML"` uses
  conditional sum of squares for starting values followed by maximum
  likelihood. `"CSS"` uses conditional sum of squares only and may be
  faster, especially for longer forecast horizons, but can produce
  different estimates. Defaults to `"CSS-ML"`.

- ...:

  Not used.

## Value

A list of class `c("acmcp", "cpforecast", "cvforecast", "forecast")`
with the following components:

- x:

  The original time series.

- series:

  The name of the series `x`.

- xreg:

  Exogenous predictor variables used, if applicable.

- method:

  A character string "acmcp".

- cp_times:

  An integer vector giving the number of conformal predictions performed
  in cross-validation for each forecast horizon.

- scorecast_times:

  An integer vector giving the number of successful scorecasts for each
  forecast horizon. Returned when `scorecast = TRUE`.

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

  A list containing information about the conformal prediction model:
  the resolved arguments in `model$args`, the call and the arguments of
  the underlying cross-validation in `model$cvforecast`, the learning
  rates actually used in `model$lr_update`, and the recursion state
  needed to extend the results in `model$state` and `model$t_last`. The
  integrator and scorecaster series are returned in `model$integrator`
  and `model$scorecaster` when `integrate` and `scorecast` are `TRUE`;
  each holds a `lower` and an `upper` element.

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
AcMCP method to capture the relationship between the \\h\\-step-ahead
forecast error and the past errors.

Scorecasts are constructed recursively, so longer forecast horizons
require more history before all scorecaster inputs are available. For
horizon \\h\\, the first scorecast can be computed at cross-validation
error index `ncal` + \\h(h-1)/2\\. The number of successful scorecasts
is reported in `scorecast_times`.

## References

Wang, X., and Hyndman, R. J. (2024). "Online conformal inference for
multi-step time series forecasting", arXiv preprint arXiv:2410.13115.

## See also

[`cvforecast`](https://xqnwang.github.io/conformalForecast/reference/cvforecast.md)
to produce `object`, and
[`update.cpforecast`](https://xqnwang.github.io/conformalForecast/reference/update.cpforecast.md)
to extend the results with new observations.

Other conformal prediction methods:
[`acp()`](https://xqnwang.github.io/conformalForecast/reference/acp.md),
[`conformal()`](https://xqnwang.github.io/conformalForecast/reference/conformal.md),
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

# AcMCP setup
Tg <- 200; delta <- 0.01
Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
KI <- 2

# AcMCP with integrator and scorecaster
acmcpfc <- acmcp(fc, ncal = 50, rolling = TRUE, KI = KI, Csat = Csat)
print(acmcpfc)
#> ACMCP 
#> 
#> Call:
#>  acmcp(object = fc, ncal = 50, rolling = TRUE, Csat = Csat, KI = KI) 
#> 
#>  cp_times (the forward step included): 101 (h=1), 100 (h=2), 99 (h=3)
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201      0.6271538 -1.013545 2.673790
#> 202      0.8607034 -2.090441 4.044708
#> 203      0.4935805 -2.055964 3.526965
summary(acmcpfc)
#> ACMCP 
#> 
#> Call:
#>  acmcp(object = fc, ncal = 50, rolling = TRUE, Csat = Csat, KI = KI) 
#> 
#>  cp_times (the forward step included): 101 (h=1), 100 (h=2), 99 (h=3)
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201      0.6271538 -1.013545 2.673790
#> 202      0.8607034 -2.090441 4.044708
#> 203      0.4935805 -2.055964 3.526965
#> 
#> Cross-validation error measures:
#>       ME   MAE   MSE RMSE    MPE    MAPE  MASE RMSSE Winkler_95 MSIS_95
#> CV 0.007 0.946 1.415 1.06 -3.933 269.763 0.992 0.882      6.618   7.132
acmcpfc$scorecast_times
#> [1] 101 100  98
```
