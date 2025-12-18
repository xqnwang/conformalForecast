# Classical split conformal prediction method

Compute prediction intervals and other information by applying the
classical split conformal prediction (SCP) method.

## Usage

``` r
scp(
  object,
  alpha = 1 - 0.01 * object$level,
  symmetric = FALSE,
  ncal = 10,
  rolling = FALSE,
  quantiletype = 1,
  weightfun = NULL,
  kess = FALSE,
  update = FALSE,
  na.rm = TRUE,
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

  Length of the calibration set. If `rolling = FALSE`, it denotes the
  initial period of calibration sets. Otherwise, it indicates the period
  of every rolling calibration set.

- rolling:

  If `TRUE`, a rolling window strategy will be adopted to form the
  calibration set. Otherwise, expanding window strategy will be used.

- quantiletype:

  An integer between 1 and 9 determining the type of quantile estimator
  to be used. Types 1 to 3 are for discontinuous quantiles, types 4 to 9
  are for continuous quantiles. See the
  [`weighted_quantile`](https://mjskay.github.io/ggdist/reference/weighted_quantile.html)
  function in the ggdist package.

- weightfun:

  Function to return a vector of weights used for sample quantile
  computation. Its first argument must be an integer indicating the
  number of observations for which weights are generated. If `NULL`,
  equal weights will be used for sample quantile computation. Currently,
  only non-data-dependent weights are supported.

- kess:

  If `TRUE`, Kish's effective sample size is used for sample quantile
  computation.

- update:

  If `TRUE`, the function will be compatible with the
  `update`([update.cpforecast](https://xqnwang.github.io/conformalForecast/reference/update.cpforecast.md))
  function, allowing for easy updates of conformal prediction.

- na.rm:

  If `TRUE`, corresponding entries in sample values and weights are
  removed if either is `NA` when calculating sample quantile.

- ...:

  Other arguments are passed to `weightfun`.

## Value

A list of class `c("scp", "cvforecast", "forecast")` with the following
components:

- x:

  The original time series.

- series:

  The name of the series `x`.

- xreg:

  Exogenous predictor variables used, if applicable.

- method:

  A character string "scp".

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

  A list containing detailed information about the `cvforecast` and
  `conformal` models.

If `mean` is included in the `object`, the components `mean`, `lower`,
and `upper` will also be returned, showing the information about the
test set forecasts generated using all available observations.

## Details

Consider a vector \\s\_{t+h\|t}\\ that contains the nonconformity scores
for the \\h\\-step-ahead forecasts.

If `symmetric` is `TRUE`, \\s\_{t+h\|t}=\|e\_{t+h\|t}\|\\. When
`rolling` is `FALSE`, the \\(1-\alpha)\\-quantile \\\hat{q}\_{t+h\|t}\\
are computed successively on expanding calibration sets
\\s\_{1+h\|1},\dots,s\_{t\|t-h}\\, for \\t=\mathrm{ncal}+h,\dots,T\\.
Then the prediction intervals will be
\\\[\hat{y}\_{t+h\|t}-\hat{q}\_{t+h\|t},
\hat{y}\_{t+h\|t}+\hat{q}\_{t+h\|t}\]\\. When `rolling` is `TRUE`, the
calibration sets will be of same length `ncal`.

If `symmetric` is `FALSE`, \\s\_{t+h\|t}^{u}=e\_{t+h\|t}\\ for upper
interval bounds and \\s\_{t+h\|t}^{l} = -e\_{t+h\|t}\\ for lower bounds.
Instead of computing \\(1-\alpha)\\-quantile, \\(1-\alpha/2)\\-quantiles
for lower bound (\\\hat{q}\_{t+h\|t}^{l}\\) and upper bound
(\\\hat{q}\_{t+h\|t}^{u}\\) are calculated based on their nonconformity
scores, respectively. Then the prediction intervals will be
\\\[\hat{y}\_{t+h\|t}-\hat{q}\_{t+h\|t}^{l},
\hat{y}\_{t+h\|t}+\hat{q}\_{t+h\|t}^{u}\]\\.

## See also

[`weighted_quantile`](https://mjskay.github.io/ggdist/reference/weighted_quantile.html)

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
scpfc <- scp(fc, symmetric = FALSE, ncal = 50, rolling = TRUE)
print(scpfc)
#> SCP 
#> 
#> Call:
#>  scp(object = fc, symmetric = FALSE, ncal = 50, rolling = TRUE) 
#> 
#>  cp_times = 99 (the forward step included) 
#> 
#> Forecasts of the forward step:
#> Cross-validation
#> 
#> Call:
#>  scp(object = fc, symmetric = FALSE, ncal = 50, rolling = TRUE) 
#> 
#>  fit_times =  (the forward step included) 
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201     -0.6049760 -2.837393 2.103852
#> 202     -0.5913297 -2.376767 1.900967
#> 203     -0.1246454 -1.950841 2.334243
summary(scpfc)
#> SCP 
#> 
#> Call:
#>  scp(object = fc, symmetric = FALSE, ncal = 50, rolling = TRUE) 
#> 
#>  cp_times = 99 (the forward step included) 
#> 
#> Forecasts of the forward step:
#> Cross-validation
#> 
#> Call:
#>  scp(object = fc, symmetric = FALSE, ncal = 50, rolling = TRUE) 
#> 
#>  fit_times =  (the forward step included) 
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201     -0.6049760 -2.837393 2.103852
#> 202     -0.5913297 -2.376767 1.900967
#> 203     -0.1246454 -1.950841 2.334243
#> 
#> Cross-validation error measures:
#>       ME   MAE   MSE  RMSE    MPE    MAPE  MASE RMSSE Winkler_95 MSIS_95
#> CV 0.077 1.036 1.672 1.169 -69.16 546.408 1.112 0.973       6.45   6.492

# Classical conformal prediction with exponential weights
expweight <- function(n) {
  0.99^{n+1-(1:n)}
}
scpfc_exp <- scp(fc, symmetric = FALSE, ncal = 50, rolling = TRUE,
                 weightfun = expweight, kess = TRUE)
```
