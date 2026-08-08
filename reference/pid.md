# Conformal PID control method

Compute prediction intervals and other information by applying the
conformal PID (Proportional-Integral-Derivative) control method.

## Usage

``` r
pid(
  object,
  alpha = 1 - 0.01 * object$level,
  symmetric = FALSE,
  ncal = 10,
  rolling = FALSE,
  integrate = TRUE,
  scorecast = !symmetric,
  scorecastfun = NULL,
  lr = 0.1,
  Tg = NROW(object$ERROR),
  delta = 0.01,
  Csat = NULL,
  KI = max(abs(object$ERROR), na.rm = TRUE),
  update = FALSE,
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

- symmetric:

  If `TRUE`, symmetric nonconformity scores (i.e. \\\|e\_{t+h\|t}\|\\)
  are used. If `FALSE`, asymmetric nonconformity scores (i.e.
  \\e\_{t+h\|t}\\) are used, and then upper bounds and lower bounds are
  produced separately. Defaults to `FALSE`.

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

  If `TRUE`, scorecasting will be included in the update process, and
  `scorecastfun` should be given. Defaults to `!symmetric`.

- scorecastfun:

  A scorecaster function to return an object of class `forecast`. Its
  first argument must be a univariate time series, and it must have an
  argument `h` for the forecast horizon. Defaults to `NULL`, which is
  only allowed when `scorecast = FALSE`.

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

- ...:

  Other arguments are passed to the `scorecastfun` function.

## Value

A list of class `c("pid", "cpforecast", "cvforecast", "forecast")` with
the following components:

- x:

  The original time series.

- series:

  The name of the series `x`.

- xreg:

  Exogenous predictor variables used, if applicable.

- method:

  A character string "pid".

- cp_times:

  An integer vector giving the number of conformal predictions performed
  in cross-validation for each forecast horizon.

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
  when `symmetric = FALSE` each of the two holds a `lower` and an
  `upper` element.

If `mean` is included in the `object`, the components `mean`, `lower`,
and `upper` will also be returned, showing the information about the
forecasts generated using all available observations.

## Details

The PID method combines three modules to make the final iteration:
\$\$q\_{t+h\|t}=\underbrace{q\_{t+h-1\|t-1} +
\eta(\mathrm{err}\_{t\|t-h}-\alpha)}\_{\mathrm{P}}+\underbrace{r_t\left(\sum\_{i=1}^t\left(\mathrm{err}\_{i\|i-h}-\alpha\right)\right)}\_{\mathrm{I}}+\underbrace{\hat{s}\_{t+h\|t}}\_{\mathrm{D}}\$\$
for each individual forecast horizon `h`, respectively, where

- Quantile tracking part (P) is \\q\_{t+h-1\|t-1} +
  \eta(\mathrm{err}\_{t\|t-h}-\alpha)\\, where \\q\_{1+h\|1}\\ is set to
  0 without a loss of generality, \\\mathrm{err}\_{t\|t-h}=1\\ if
  \\s\_{t\|t-h}\>q\_{t\|t-h}\\, and \\\mathrm{err}\_{t\|t-h}=0\\ if
  \\s\_{t\|t-h} \leq q\_{t\|t-h}\\.

- Error integration part (I) is
  \\r_t\left(\sum\_{i=1}^t\left(\mathrm{err}\_{i\|i-h}-\alpha\right)\right)\\.
  Here we use a nonlinear saturation function \\r_t(x)=K\_{\mathrm{I}}
  \tan \left(x \log (t) /\left(t C\_{\text {sat }}\right)\right)\\,
  where we set \\\tan (x)=\operatorname{sign}(x) \cdot \infty\\ for \\x
  \notin\[-\pi / 2, \pi / 2\]\\, and \\C\_{\text {sat }},
  K\_{\mathrm{I}}\>0\\ are constants that we choose heuristically.

- Scorecasting part (D) is \\\hat{s}\_{t+h\|t}\\, the forecast generated
  by training a scorecaster based on the nonconformity scores available
  at time \\t\\.

## References

Angelopoulos, A., Candes, E., and Tibshirani, R. J. (2023). "Conformal
PID control for time series prediction", *Advances in Neural Information
Processing Systems*, **36**, 23047–23074.

## See also

[`cvforecast`](https://xqnwang.github.io/conformalForecast/reference/cvforecast.md)
to produce `object`, and
[`update.cpforecast`](https://xqnwang.github.io/conformalForecast/reference/update.cpforecast.md)
to extend the results with new observations.

Other conformal prediction methods:
[`acmcp()`](https://xqnwang.github.io/conformalForecast/reference/acmcp.md),
[`acp()`](https://xqnwang.github.io/conformalForecast/reference/acp.md),
[`conformal()`](https://xqnwang.github.io/conformalForecast/reference/conformal.md),
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

# PID setup
Tg <- 200; delta <- 0.01
Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
KI <- 2

# PID without scorecaster
pidfc_nsf <- pid(fc, ncal = 50, rolling = TRUE, scorecast = FALSE,
                 KI = KI, Csat = Csat)
print(pidfc_nsf)
#> PID 
#> 
#> Call:
#>  pid(object = fc, ncal = 50, rolling = TRUE, scorecast = FALSE,  
#>      Csat = Csat, KI = KI) 
#> 
#>  cp_times (the forward step included): 101 (h=1), 100 (h=2), 99 (h=3)
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201      0.6271538 -1.189121 3.134690
#> 202      0.8607034 -1.600259 3.730935
#> 203      0.4935805 -2.173248 3.218838
summary(pidfc_nsf)
#> PID 
#> 
#> Call:
#>  pid(object = fc, ncal = 50, rolling = TRUE, scorecast = FALSE,  
#>      Csat = Csat, KI = KI) 
#> 
#>  cp_times (the forward step included): 101 (h=1), 100 (h=2), 99 (h=3)
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201      0.6271538 -1.189121 3.134690
#> 202      0.8607034 -1.600259 3.730935
#> 203      0.4935805 -2.173248 3.218838
#> 
#> Cross-validation error measures:
#>       ME   MAE   MSE RMSE    MPE    MAPE  MASE RMSSE Winkler_95 MSIS_95
#> CV 0.007 0.946 1.415 1.06 -3.933 269.763 0.992 0.882      6.339    6.82

# PID with a Naive model for the scorecaster
naivefun <- function(x, h) {
  naive(x) |> forecast(h = h)
}
pidfc <- pid(fc, ncal = 50, rolling = TRUE, scorecastfun = naivefun,
             KI = KI, Csat = Csat)
```
