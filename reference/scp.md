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

  Length of the calibration set. If `rolling = FALSE`, it denotes the
  initial period of the calibration sets. Otherwise, it indicates the
  period of every rolling calibration set. Defaults to `10`.

- rolling:

  If `TRUE`, a rolling window strategy will be adopted to form the
  calibration set. Otherwise, an expanding window strategy will be used.
  Defaults to `FALSE`.

- quantiletype:

  An integer between 1 and 9 determining the type of quantile estimator
  to be used. Types 1 to 3 are for discontinuous quantiles, types 4 to 9
  are for continuous quantiles. See the
  [`weighted_quantile`](https://mjskay.github.io/ggdist/reference/weighted_quantile.html)
  function in the ggdist package. Defaults to `1`.

- weightfun:

  Function to return a vector of weights used for the sample quantile
  computation. Its first argument must be an integer indicating the
  number of observations for which weights are generated. Defaults to
  `NULL`, in which case equal weights are used. Currently, only
  non-data-dependent weights are supported.

- kess:

  If `TRUE`, Kish's effective sample size is used for the sample
  quantile computation. Defaults to `FALSE`.

- update:

  If `TRUE`, `object` already holds the results of a previous call and
  only the newly added time steps are computed; the prediction intervals
  produced earlier are carried over unchanged. Set by
  [`update.cpforecast`](https://xqnwang.github.io/conformalForecast/reference/update.cpforecast.md)
  and not normally set by hand. Defaults to `FALSE`.

- na.rm:

  If `TRUE`, corresponding entries in the sample values and weights are
  removed if either is `NA` when calculating the sample quantile.
  Defaults to `TRUE`.

- ...:

  Other arguments are passed to `weightfun`.

## Value

A list of class `c("scp", "cpforecast", "cvforecast", "forecast")` with
the following components:

- x:

  The original time series.

- series:

  The name of the series `x`.

- xreg:

  Exogenous predictor variables used, if applicable.

- method:

  A character string "scp".

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
  the resolved arguments in `model$args`, and the call and the arguments
  of the underlying cross-validation in `model$cvforecast`.

If `mean` is included in the `object`, the components `mean`, `lower`,
and `upper` will also be returned, showing the information about the
test set forecasts generated using all available observations.

## Details

Consider a vector \\s\_{t+h\|t}\\ that contains the nonconformity scores
for the \\h\\-step-ahead forecasts.

If `symmetric` is `TRUE`, \\s\_{t+h\|t}=\|e\_{t+h\|t}\|\\. When
`rolling` is `FALSE`, the \\(1-\alpha)\\-quantiles \\\hat{q}\_{t+h\|t}\\
are computed successively on expanding calibration sets
\\s\_{1+h\|1},\dots,s\_{t\|t-h}\\, for \\t=\mathrm{ncal}+h,\dots,T\\.
Then the prediction intervals will be
\\\[\hat{y}\_{t+h\|t}-\hat{q}\_{t+h\|t},
\hat{y}\_{t+h\|t}+\hat{q}\_{t+h\|t}\]\\. When `rolling` is `TRUE`, the
calibration sets will be of the same length `ncal`.

If `symmetric` is `FALSE`, \\s\_{t+h\|t}^{u}=e\_{t+h\|t}\\ for the upper
interval bounds and \\s\_{t+h\|t}^{l} = -e\_{t+h\|t}\\ for the lower
bounds. Instead of computing the \\(1-\alpha)\\-quantile, the
\\(1-\alpha/2)\\-quantiles for the lower bound \\\hat{q}\_{t+h\|t}^{l}\\
and the upper bound \\\hat{q}\_{t+h\|t}^{u}\\ are calculated based on
their nonconformity scores, respectively. Then the prediction intervals
will be \\\[\hat{y}\_{t+h\|t}-\hat{q}\_{t+h\|t}^{l},
\hat{y}\_{t+h\|t}+\hat{q}\_{t+h\|t}^{u}\]\\.

## See also

[`cvforecast`](https://xqnwang.github.io/conformalForecast/reference/cvforecast.md)
to produce `object`,
[`update.cpforecast`](https://xqnwang.github.io/conformalForecast/reference/update.cpforecast.md)
to extend the results with new observations, and
[`weighted_quantile`](https://mjskay.github.io/ggdist/reference/weighted_quantile.html)
for the weighted quantile estimator used here.

Other conformal prediction methods:
[`acmcp()`](https://xqnwang.github.io/conformalForecast/reference/acmcp.md),
[`acp()`](https://xqnwang.github.io/conformalForecast/reference/acp.md),
[`conformal()`](https://xqnwang.github.io/conformalForecast/reference/conformal.md),
[`pid()`](https://xqnwang.github.io/conformalForecast/reference/pid.md)

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
scpfc <- scp(fc, ncal = 50, rolling = TRUE)
print(scpfc)
#> SCP 
#> 
#> Call:
#>  scp(object = fc, ncal = 50, rolling = TRUE) 
#> 
#>  cp_times (the forward step included): 101 (h=1), 100 (h=2), 99 (h=3)
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201      0.6271538 -1.023525 3.234253
#> 202      0.8607034 -1.291456 4.064625
#> 203      0.4935805 -1.675089 3.691253
summary(scpfc)
#> SCP 
#> 
#> Call:
#>  scp(object = fc, ncal = 50, rolling = TRUE) 
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

# Classical conformal prediction with exponential weights
expweight <- function(n) {
  0.99^{n+1-(1:n)}
}
scpfc_exp <- scp(fc, ncal = 50, rolling = TRUE,
                 weightfun = expweight, kess = TRUE)
```
