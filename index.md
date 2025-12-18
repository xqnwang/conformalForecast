# conformalForecast

The R package *conformalForecast* provides methods and tools for
performing multistep-ahead time series forecasting using conformal
prediction methods including classical conformal prediction, adaptive
conformal prediction, conformal PID control, and autocorrelated
multistep-ahead conformal prediction.

## Installation

You can install the development version of conformalForecast from
[GitHub](https://github.com/xqnwang/conformalForecast) with:

``` r
# install.packages("remotes")
remotes::install_github("xqnwang/conformalForecast")
```

You can also get the official release version from CRAN:

``` r
install.packages("conformalForecast")
```

## Example

This is a basic example which shows you how to perform a classical
conformal prediction method:

``` r
library(conformalForecast)
library(forecast)

# Simulate time series from an AR(2) model
series <- arima.sim(n = 1000, list(ar = c(0.8, -0.5)), sd = sqrt(1))

# Time series cross-validation
far2 <- function(x, h, level) {
  Arima(x, order = c(2, 0, 0)) |>
    forecast(h = h, level)
}
fc <- cvforecast(series, forecastfun = far2, h = 3, level = c(80, 95),
                 forward = TRUE, initial = 1, window = 100)

# Classical conformal prediction
scpfc <- scp(fc, symmetric = FALSE, ncal = 100, rolling = TRUE)

# Interval forecast accuracy
accuracy(scpfc, byhorizon = TRUE)

# Mean coverage
coverage(scpfc, window = 500, level = 95)

# Mean and median interval width
width(scpfc, window = 500, level = 95, includemedian = TRUE)
```

## License

This package is free and open source software, licensed under GPL-3.
