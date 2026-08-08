# Calculate interval forecast coverage

Calculate the mean coverage and the `ifinn` matrix for the prediction
intervals on the validation set. If `window` is not `NULL`, a matrix of
the rolling means of the interval forecast coverage is also returned.

## Usage

``` r
coverage(object, ..., level = 95, window = NULL, na.rm = FALSE)
```

## Arguments

- object:

  An object of class `"cvforecast"` or `"cpforecast"`.

- ...:

  Arguments `x`, `LOWER`, and `UPPER` if `object` is missing. Bounds may
  be time-series matrices or lists keyed by confidence level.

- level:

  Target confidence level for the prediction intervals. Only one level
  can be specified. Defaults to `95`.

- window:

  If not `NULL`, the rolling mean matrix for the coverage is also
  returned. Defaults to `NULL`.

- na.rm:

  A logical indicating whether `NA` values should be stripped before the
  rolling mean computation proceeds. Defaults to `FALSE`.

## Value

A list of class `"coverage"` with the following components:

- mean:

  Mean coverage across the validation set.

- ifinn:

  An indicator matrix as a multivariate time series, where the \\h\\th
  column holds the coverage for forecast horizon \\h\\. The time index
  corresponds to the period for which the forecast is produced.

- rollmean:

  If `window` is not `NULL`, a matrix of the rolling means of the
  interval forecast coverage will be returned.

## See also

[`cvforecast`](https://xqnwang.github.io/conformalForecast/reference/cvforecast.md)
and the conformal methods, which produce `object`.

Other evaluation functions:
[`accuracy.cvforecast()`](https://xqnwang.github.io/conformalForecast/reference/accuracy.cvforecast.md),
[`width()`](https://xqnwang.github.io/conformalForecast/reference/width.md)

## Examples

``` r
# Simulate time series from an AR(2) model
library(forecast)
set.seed(1)
series <- arima.sim(n = 200, list(ar = c(0.8, -0.5)), sd = sqrt(1))

# Cross-validation forecasting with a rolling window
far2 <- function(x, h, level) {
  Arima(x, order = c(2, 0, 0)) |>
    forecast(h = h, level)
}
fc <- cvforecast(series, forecastfun = far2, h = 3, level = 95, window = 50)

# Mean and rolling mean coverage for interval forecasts on validation set
cov_fc <- coverage(fc, level = 95, window = 50)
str(cov_fc)
#> List of 3
#>  $ mean    : Named num [1:3] 0.927 0.906 0.912
#>   ..- attr(*, "names")= chr [1:3] "h=1" "h=2" "h=3"
#>  $ ifinn   : Time-Series [1:150, 1:3] from 51 to 200: TRUE TRUE TRUE TRUE TRUE TRUE ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "h=1" "h=2" "h=3"
#>  $ rollmean: Time-Series [1:101, 1:3] from 100 to 200: 0.98 0.98 0.98 0.98 0.98 0.98 0.98 0.98 0.98 0.98 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "h=1" "h=2" "h=3"
#>  - attr(*, "class")= chr "coverage"

# Coverage calculated directly from interval components
coverage(x = fc$x, LOWER = fc$LOWER, UPPER = fc$UPPER, level = 95)
#>       h=1       h=2       h=3 
#> 0.9266667 0.9060403 0.9121622 
```
