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

- symmetric:

  If `TRUE`, symmetric nonconformity scores (i.e. \\\|e\_{t+h\|t}\|\\)
  are used. If `FALSE`, asymmetric nonconformity scores (i.e.
  \\e\_{t+h\|t}\\) are used, and then upper bounds and lower bounds are
  produced separately.

- ncal:

  Length of the burn-in period for training the scorecaster. If
  `rolling = TRUE`, it is also used as the length of the trailing
  windows for learning rate calculation and the windows for the
  calibration set. If `rolling = FALSE`, it is used as initial period of
  calibration sets and trailing windows for learning rate calculation.

- rolling:

  If `TRUE`, a rolling window strategy will be adopted to form the
  trailing window for learning rate calculation and the calibration set
  for scorecaster if applicable. Otherwise, expanding window strategy
  will be used.

- integrate:

  If `TRUE`, error integration will be included in the update process.

- scorecast:

  If `TRUE`, scorecasting will be included in the update process, and
  `scorecastfun` should be given.

- scorecastfun:

  A scorecaster function to return an object of class `forecast`. Its
  first argument must be a univariate time series, and it must have an
  argument `h` for the forecast horizon.

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

  Other arguments are passed to the `scorecastfun` function.

## Value

A list of class `c("pid", "cpforecast", "forecast")` with the following
components:

- x:

  The original time series.

- series:

  The name of the series `x`.

- method:

  A character string "pid".

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

- Scorecasting part (D) is \\\hat{s}\_{t+h\|t}\\ is forecast generated
  by training a scorecaster based on nonconformity scores available at
  time \\t\\.

## References

Angelopoulos, A., Candes, E., and Tibshirani, R. J. (2024). "Conformal
PID control for time series prediction", *Advances in Neural Information
Processing Systems*, **36**, 23047–23074.

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
# PID setup
Tg <- 200; delta <- 0.01
Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
KI <- 2
lr <- 0.1
# PID without scorecaster
pidfc_nsf <- pid(fc, symmetric = FALSE, ncal = 50, rolling = TRUE,
                 integrate = TRUE, scorecast = FALSE,
                 lr = lr, KI = KI, Csat = Csat)
print(pidfc_nsf)
#> PID 
#> 
#> Call:
#>  pid(object = fc, symmetric = FALSE, ncal = 50, rolling = TRUE,  
#>      integrate = TRUE, scorecast = FALSE, lr = lr, Csat = Csat,  
#>      KI = KI) 
#> 
#>  cp_times = 148 (the forward step included) 
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201     0.78863033 -1.494426 2.585959
#> 202     0.19240047 -1.821271 2.751143
#> 203    -0.02560772 -1.919723 2.370350
summary(pidfc_nsf)
#> PID 
#> 
#> Call:
#>  pid(object = fc, symmetric = FALSE, ncal = 50, rolling = TRUE,  
#>      integrate = TRUE, scorecast = FALSE, lr = lr, Csat = Csat,  
#>      KI = KI) 
#> 
#>  cp_times = 148 (the forward step included) 
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201     0.78863033 -1.494426 2.585959
#> 202     0.19240047 -1.821271 2.751143
#> 203    -0.02560772 -1.919723 2.370350
#> 
#> Cross-validation error measures:
#>       ME   MAE   MSE  RMSE     MPE    MAPE  MASE RMSSE Winkler_95 MSIS_95
#> CV 0.144 0.979 1.508 1.105 -89.391 522.845 0.993 0.866      6.066   5.969
# PID with a Naive model for the scorecaster
naivefun <- function(x, h) {
  naive(x) |> forecast(h = h)
}
pidfc <- pid(fc, symmetric = FALSE, ncal = 50, rolling = TRUE,
             integrate = TRUE, scorecast = TRUE, scorecastfun = naivefun,
             lr = lr, KI = KI, Csat = Csat)
```
