
<!-- README.md is generated from README.Rmd. Please edit that file -->

# conformalForecast

<!-- badges: start -->

[![R-CMD-check](https://github.com/xqnwang/conformalForecast/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/xqnwang/conformalForecast/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/conformalForecast)](https://CRAN.R-project.org/package=conformalForecast)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/grand-total/conformalForecast)](https://cran.r-project.org/package=conformalForecast)
[![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://github.com/xqnwang/conformalForecast/blob/main/LICENSE.md)
<!-- badges: end -->

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

For a detailed tutorial and examples, please see the
[vignette](https://xqnwang.github.io/conformalForecast/articles/conformalForecast.html).

This basic example shows how to perform classical conformal prediction
and update the results when new observations become available:

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
scpfc <- conformal(fc, method = "scp", symmetric = FALSE, ncal = 100,
                   rolling = TRUE)

# Update conformal prediction with newly available observations
scpfc_updated <- update(scpfc, new_data = c(1.5, 0.8, 2.3),
                        forecastfun = far2)

# Interval forecast accuracy
accuracy(scpfc_updated, byhorizon = TRUE)

# Mean coverage
coverage(scpfc_updated, window = 500, level = 95)

# Mean and median interval width
width(scpfc_updated, window = 500, level = 95, includemedian = TRUE)
```

## License

This package is free and open source software, licensed under GPL-3.
