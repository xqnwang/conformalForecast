# Time series cross-validation forecasting

Compute forecasts and other information by applying `forecastfun` to
subsets of the time series `y` using a rolling forecast origin.

## Usage

``` r
cvforecast(
  y,
  forecastfun,
  h = 1,
  level = c(80, 95),
  forward = TRUE,
  xreg = NULL,
  initial = 1,
  window = NULL,
  ...
)
```

## Arguments

- y:

  Univariate time series.

- forecastfun:

  Function to return an object of class `"forecast"`. Its first argument
  must be a univariate time series, and it must have an argument `h` for
  the forecast horizon and an argument `level` for the confidence level
  for prediction intervals. If exogenous predictors are used, then it
  must also have `xreg` and `newxreg` arguments corresponding to the
  training and test periods, respectively.

- h:

  Forecast horizon.

- level:

  Confidence level for prediction intervals. If `NULL`, prediction
  intervals will not be generated.

- forward:

  If `TRUE`, the final forecast origin for forecasting is \\y_T\\.
  Otherwise, the final forecast origin is \\y\_{T-1}\\.

- xreg:

  Exogenous predictor variables passed to `forecastfun` if required. It
  should be of the same size as `y`+`forward`\*`h`, otherwise, `NA`
  padding or subsetting will be applied.

- initial:

  Initial period of the time series where no cross-validation
  forecasting is performed.

- window:

  Length of the rolling window. If `NULL`, a rolling window will not be
  used.

- ...:

  Other arguments are passed to `forecastfun`.

## Value

A list of class `c("cvforecast", "forecast")` with components:

- x:

  The original time series.

- series:

  The name of the series `x`.

- xreg:

  Exogenous predictor variables used in the model, if applicable.

- method:

  A character string "cvforecast".

- fit_times:

  The number of times the model is fitted in cross-validation.

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

- forward:

  Whether `forward` is applied.

If `forward` is `TRUE`, the components `mean`, `lower`, `upper`, and
`model` will also be returned, showing the information about the final
fitted model and forecasts using all available observations, see e.g.
[`forecast.ets`](https://pkg.robjhyndman.com/forecast/reference/forecast.ets.html)
for more details.

## Details

Let `y` denote the time series \\y_1,\dots,y_T\\ and let \\t_0\\ denote
the initial period.

Suppose `forward = TRUE`. If `window` is `NULL`, `forecastfun` is
applied successively to the subset time series \\y\_{1},\dots,y_t\\, for
\\t=t_0,\dots,T\\, generating forecasts
\\\hat{y}\_{t+1\|t},\dots,\hat{y}\_{t+h\|t}\\. If `window` is not `NULL`
and has a length of \\t_w\\, then `forecastfun` is applied successively
to the subset time series \\y\_{t-t_w+1},\dots,y\_{t}\\, for
\\t=\max(t_0, t_w),\dots,T\\.

If `forward` is `FALSE`, the last observation used for training will be
\\y\_{T-1}\\.

## Examples

``` r
# Simulate time series from an AR(2) model
library(forecast)
series <- arima.sim(n = 200, list(ar = c(0.8, -0.5)), sd = sqrt(1))

# Example with a rolling window
far2 <- function(x, h, level) {
  Arima(x, order = c(2, 0, 0)) |>
    forecast(h = h, level)
}
fc <- cvforecast(series, forecastfun = far2, h = 3, level = 95,
                 forward = TRUE, initial = 1, window = 50)
print(fc)
#> Cross-validation
#> 
#> Call:
#>  cvforecast(y = series, forecastfun = far2, h = 3, level = 95,  
#>      forward = TRUE, initial = 1, window = 50) 
#> 
#>  fit_times = 151 (the forward step included) 
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201     -1.1554197 -2.740135 0.429296
#> 202      0.2822432 -1.902485 2.466971
#> 203      1.0323863 -1.174433 3.239206
summary(fc)
#> Cross-validation
#> 
#> Call:
#>  cvforecast(y = series, forecastfun = far2, h = 3, level = 95,  
#>      forward = TRUE, initial = 1, window = 50) 
#> 
#>  fit_times = 151 (the forward step included) 
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 201     -1.1554197 -2.740135 0.429296
#> 202      0.2822432 -1.902485 2.466971
#> 203      1.0323863 -1.174433 3.239206
#> 
#> Cross-validation error measures:
#>       ME   MAE   MSE  RMSE     MPE    MAPE  MASE RMSSE Winkler_95 MSIS_95
#> CV -0.02 0.997 1.552 1.112 -18.514 291.627 0.941  0.84      5.526     5.2

# Example with exogenous predictors
far2_xreg <- function(x, h, level, xreg, newxreg) {
  Arima(x, order=c(2, 0, 0), xreg = xreg) |>
    forecast(h = h, level = level, xreg = newxreg)
}
fc_xreg <- cvforecast(series, forecastfun = far2_xreg, h = 3, level = 95,
                      forward = TRUE, xreg = matrix(rnorm(406), ncol = 2, nrow = 203),
                      initial = 1, window = 50)
```
