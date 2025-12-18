# Adaptive conformal prediction method

Compute prediction intervals and other information by applying the
adaptive conformal prediction (ACP) method.

## Usage

``` r
acp(
  object,
  alpha = 1 - 0.01 * object$level,
  gamma = 0.005,
  symmetric = FALSE,
  ncal = 10,
  rolling = FALSE,
  quantiletype = 1,
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

- gamma:

  The step size parameter \\\gamma\>0\\ for \\\alpha\\ updating.

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

- update:

  If `TRUE`, the function will be compatible with the
  `update`([update.cpforecast](https://xqnwang.github.io/conformalForecast/reference/update.cpforecast.md))
  function, allowing for easy updates of conformal prediction.

- na.rm:

  If `TRUE`, corresponding entries in sample values are removed if it is
  `NA` when calculating sample quantile.

- ...:

  Other arguments are passed to the
  [`weighted_quantile`](https://mjskay.github.io/ggdist/reference/weighted_quantile.html)
  function for quantile computation.

## Value

A list of class `c("acp", "cpforecast", "forecast")` with the following
components:

- x:

  The original time series.

- series:

  The name of the series `x`.

- method:

  A character string "acp".

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

The ACP method considers the online update:
\$\$\alpha\_{t+h\|t}:=\alpha\_{t+h-1\|t-1}+\gamma(\alpha-\mathrm{err}\_{t\|t-h}),\$\$
for each individual forecast horizon `h`, respectively, where
\\\mathrm{err}\_{t\|t-h}=1\\ if \\s\_{t\|t-h}\>q\_{t\|t-h}\\, and
\\\mathrm{err}\_{t\|t-h}=0\\ if \\s\_{t\|t-h} \leq q\_{t\|t-h}\\.

## References

Gibbs, I., and Candes, E. (2021). "Adaptive conformal inference under
distribution shift", *Advances in Neural Information Processing
Systems*, **34**, 1660–1672.

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

# ACP with asymmetric nonconformity scores and rolling calibration sets
acpfc <- acp(fc, symmetric = FALSE, gamma = 0.005, ncal = 50, rolling = TRUE)
print(acpfc)
#> ACP 
#> 
#> Call:
#>  acp(object = fc, gamma = 0.005, symmetric = FALSE, ncal = 50,  
#>      rolling = TRUE) 
#> 
#>  cp_times = 99 (the forward step included) 
#> 
#> Forecasts of the forward step:
#> Cross-validation
#> 
#> Call:
#>  acp(object = fc, gamma = 0.005, symmetric = FALSE, ncal = 50,  
#>      rolling = TRUE) 
#> 
#>  fit_times =  (the forward step included) 
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201      0.1969684 -1.055741 1.795689
#> 202     -0.5275402 -2.526663 1.438041
#> 203     -0.6180594 -2.402479 1.119214
summary(acpfc)
#> ACP 
#> 
#> Call:
#>  acp(object = fc, gamma = 0.005, symmetric = FALSE, ncal = 50,  
#>      rolling = TRUE) 
#> 
#>  cp_times = 99 (the forward step included) 
#> 
#> Forecasts of the forward step:
#> Cross-validation
#> 
#> Call:
#>  acp(object = fc, gamma = 0.005, symmetric = FALSE, ncal = 50,  
#>      rolling = TRUE) 
#> 
#>  fit_times =  (the forward step included) 
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201      0.1969684 -1.055741 1.795689
#> 202     -0.5275402 -2.526663 1.438041
#> 203     -0.6180594 -2.402479 1.119214
#> 
#> Cross-validation error measures:
#>        ME  MAE   MSE RMSE      MPE    MAPE  MASE RMSSE Winkler_95 MSIS_95
#> CV -0.038 0.93 1.341 1.04 -220.016 425.573 0.958 0.875      6.023   6.002
```
