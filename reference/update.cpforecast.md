# Update and repeform cross-validation forecasting and conformal prediction

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

  A vector of newly available data.

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
  argument in `object`.

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
#>  scp(object = object, symmetric = FALSE, ncal = 50, rolling = TRUE,  
#>      update = TRUE) 
#> 
#>  cp_times = 102 (the forward step included) 
#> 
#> Forecasts of the forward step:
#> Cross-validation
#> 
#> Call:
#>  scp(object = object, symmetric = FALSE, ncal = 50, rolling = TRUE,  
#>      update = TRUE) 
#> 
#>  fit_times =  (the forward step included) 
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 204      1.2910976 -1.249629 5.592616
#> 205     -0.3496234 -3.679947 4.522049
#> 206     -0.9708415 -4.416003 3.579224
summary(scpfc_update)
#> SCP 
#> 
#> Call:
#>  scp(object = object, symmetric = FALSE, ncal = 50, rolling = TRUE,  
#>      update = TRUE) 
#> 
#>  cp_times = 102 (the forward step included) 
#> 
#> Forecasts of the forward step:
#> Cross-validation
#> 
#> Call:
#>  scp(object = object, symmetric = FALSE, ncal = 50, rolling = TRUE,  
#>      update = TRUE) 
#> 
#>  fit_times =  (the forward step included) 
#> 
#> Forecasts of the forward step:
#>     Point Forecast     Lo 95    Hi 95
#> 204      1.2910976 -1.249629 5.592616
#> 205     -0.3496234 -3.679947 4.522049
#> 206     -0.9708415 -4.416003 3.579224
#> 
#> Cross-validation error measures:
#>        ME   MAE   MSE  RMSE    MPE    MAPE MASE RMSSE Winkler_95 MSIS_95
#> CV -0.001 1.131 2.037 1.276 26.357 239.353 1.04 0.943      8.998   7.994
```
