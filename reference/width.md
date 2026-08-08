# Calculate interval forecast width

Calculate the mean width of the prediction intervals on the validation
set. If `window` is not `NULL`, a matrix of the rolling means of the
interval width is also returned. If `includemedian` is `TRUE`, the
information of the median interval width will be returned.

## Usage

``` r
width(
  object,
  ...,
  level = 95,
  includemedian = FALSE,
  window = NULL,
  na.rm = FALSE
)
```

## Arguments

- object:

  An object of class `"cvforecast"` or `"cpforecast"`.

- ...:

  Time-series matrices `LOWER` and `UPPER` if `object` is missing. They
  may also be lists keyed by confidence level.

- level:

  Target confidence level for the prediction intervals. Only one level
  can be specified. Defaults to `95`.

- includemedian:

  If `TRUE`, the median interval width will also be returned. Defaults
  to `FALSE`.

- window:

  If not `NULL`, the rolling mean (and the rolling median if applicable)
  matrix for the interval width will also be returned. Defaults to
  `NULL`.

- na.rm:

  A logical indicating whether `NA` values should be stripped before the
  rolling mean and rolling median computation proceeds. Defaults to
  `FALSE`.

## Value

A list of class `"width"` with the following components:

- width:

  Forecast interval width as a multivariate time series, where the
  \\h\\th column holds the interval width for the forecast horizon
  \\h\\. The time index corresponds to the period for which the forecast
  is produced.

- mean:

  Mean interval width across the validation set.

- rollmean:

  If `window` is not `NULL`, a matrix of the rolling means of the
  interval width will be returned.

- median:

  Median interval width across the validation set.

- rollmedian:

  If `window` is not `NULL`, a matrix of the rolling medians of the
  interval width will be returned.

## See also

[`cvforecast`](https://xqnwang.github.io/conformalForecast/reference/cvforecast.md)
and the conformal methods, which produce `object`.

Other evaluation functions:
[`accuracy.cvforecast()`](https://xqnwang.github.io/conformalForecast/reference/accuracy.cvforecast.md),
[`coverage()`](https://xqnwang.github.io/conformalForecast/reference/coverage.md)

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

# Mean and rolling mean width for interval forecasts on validation set
wid_fc <- width(fc, level = 95, window = 50)
str(wid_fc)
#> List of 3
#>  $ width   : Time-Series [1:153, 1:3] from 51 to 203: 3.66 3.67 3.64 3.6 3.59 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "h=1" "h=2" "h=3"
#>  $ mean    : Named num [1:3] 3.72 4.69 4.76
#>   ..- attr(*, "names")= chr [1:3] "h=1" "h=2" "h=3"
#>  $ rollmean: Time-Series [1:104, 1:3] from 100 to 203: 3.65 3.65 3.63 3.63 3.62 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "h=1" "h=2" "h=3"
#>  - attr(*, "class")= chr "width"

# Width calculated directly from interval components
width(LOWER = fc$LOWER, UPPER = fc$UPPER, level = 95)
#> Mean width:
#>      h=1      h=2      h=3 
#> 3.720038 4.692975 4.757742 
```
