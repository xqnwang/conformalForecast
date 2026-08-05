# Update and reperform cross-validation forecasting and conformal prediction

Update conformal prediction intervals and other information by applying
the `cvforecast` and `conformal` functions.

## Usage

``` r
# S3 method for class 'cpforecast'
update(object, new_data, forecastfun, new_xreg = NULL, ...)
```

## Arguments

- object:

  An object of class `"cpforecast"`.

- new_data:

  A non-empty numeric vector of newly available data.

- forecastfun:

  Function to return an object of class `"forecast"`. Its first argument
  must be a univariate time series, and it must have an argument `h` for
  the forecast horizon and an argument `level` for the confidence level
  for prediction intervals. If exogenous predictors are used, then it
  must also have `xreg` and `newxreg` arguments corresponding to the
  training and test periods, respectively.

- new_xreg:

  Newly available exogenous predictor variables passed to `forecastfun`
  if required. The number of rows should match the length of `new_data`,
  and the number of columns should match the dimensions of the `xreg`
  argument in `object`. These rows extend `object$xreg` and correspond
  to the periods immediately following it.

- ...:

  Other arguments are passed to `forecastfun`.

## Value

A refreshed object of class `"cpforecast"` with updated fields (e.g.,
`x`, `MEAN`, `ERROR`, `LOWER`, `UPPER`, and any method-specific
components), reflecting newly appended data and re-computed
cross-validation forecasts and conformal prediction intervals.

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
fc <- cvforecast(series, forecastfun = far2, h = 3, level = 95,
                 forward = TRUE, initial = 1, window = 50)

# Classical conformal prediction with equal weights
scpfc <- conformal(fc, method = "scp", symmetric = FALSE, ncal = 50, rolling = TRUE)

# Update conformal prediction using newly available data
scpfc_update <- update(scpfc, forecastfun = far2, new_data = c(1.5, 0.8, 2.3))
print(scpfc_update)
#> SCP 
#> 
#> Call:
#>  scp(object = object, alpha = 0.0499999999999999, symmetric = FALSE,  
#>      ncal = 50, rolling = TRUE, quantiletype = 1, weightfun = ..6,  
#>      kess = FALSE, update = TRUE, na.rm = TRUE) 
#> 
#>  cp_times (the forward step included): 104 (h=1), 103 (h=2), 102 (h=3)
#> 
#> Forecasts of the forward step:
#>     Point Forecast      Lo 95    Hi 95
#> 204      1.3171165 -0.3335623 3.753786
#> 205     -0.1244955 -2.2766551 3.079426
#> 206     -0.5533392 -2.7220088 2.644333
summary(scpfc_update)
#> SCP 
#> 
#> Call:
#>  scp(object = object, alpha = 0.0499999999999999, symmetric = FALSE,  
#>      ncal = 50, rolling = TRUE, quantiletype = 1, weightfun = ..6,  
#>      kess = FALSE, update = TRUE, na.rm = TRUE) 
#> 
#>  cp_times (the forward step included): 104 (h=1), 103 (h=2), 102 (h=3)
#> 
#> Forecasts of the forward step:
#>     Point Forecast      Lo 95    Hi 95
#> 204      1.3171165 -0.3335623 3.753786
#> 205     -0.1244955 -2.2766551 3.079426
#> 206     -0.5533392 -2.7220088 2.644333
#> 
#> Cross-validation error measures:
#>       ME   MAE   MSE  RMSE    MPE    MAPE  MASE RMSSE Winkler_95 MSIS_95
#> CV 0.037 0.957 1.448 1.072 -3.518 265.271 1.002  0.89      6.104   6.536
```
