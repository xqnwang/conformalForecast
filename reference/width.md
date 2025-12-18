# Calculate interval forecast width

Calculate the mean width of prediction intervals on the validation set.
If `window` is not `NULL`, a matrix of the rolling means of interval
width is also returned. If `includemedian` is `TRUE`, the information of
the median interval width will be returned.

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

  Additional inputs if `object` is missing.

- level:

  Target confidence level for prediction intervals.

- includemedian:

  If `TRUE`, the median interval width will also be returned.

- window:

  If not `NULL`, the rolling mean (and rolling median if applicable)
  matrix for interval width will also be returned.

- na.rm:

  A logical indicating whether `NA` values should be stripped before the
  rolling mean and rolling median computation proceeds.

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

  If `window` is not NULL, a matrix of the rolling means of interval
  width will be returned.

- median:

  Median interval width across the validation set.

- rollmedian:

  If `window` is not NULL, a matrix of the rolling medians of interval
  width will be returned.

## Examples

``` r
# Simulate time series from an AR(2) model
library(forecast)
series <- arima.sim(n = 200, list(ar = c(0.8, -0.5)), sd = sqrt(1))

# Cross-validation forecasting with a rolling window
far2 <- function(x, h, level) {
  Arima(x, order = c(2, 0, 0)) |>
    forecast(h = h, level)
}
fc <- cvforecast(series, forecastfun = far2, h = 3, level = 95,
                 forward = TRUE, initial = 1, window = 50)

# Mean and rolling mean width for interval forecasts on validation set
wid_fc <- width(fc, level = 95, window = 50)
str(wid_fc)
#> List of 3
#>  $ width   : Time-Series [1:153, 1:3] from 51 to 203: 3.87 3.96 3.94 3.95 3.95 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "h=1" "h=2" "h=3"
#>  $ mean    : Named num [1:3] 4.21 5.35 5.4
#>   ..- attr(*, "names")= chr [1:3] "h=1" "h=2" "h=3"
#>  $ rollmean: Time-Series [1:104, 1:3] from 100 to 203: 4.08 4.09 4.1 4.1 4.1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "h=1" "h=2" "h=3"
#>  - attr(*, "class")= chr "width"
```
