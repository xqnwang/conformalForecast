
# conformalForecast <img src="man/figures/logo.svg" alt="logo" align="right" width="150" style="border: none; float: right;"/>

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

You can also get the official release version from
[CRAN](https://CRAN.R-project.org/package=conformalForecast) with:

``` r
install.packages("conformalForecast")
```

## Example

This basic example shows how to perform classical conformal prediction
and update the results when new observations become available.

We simulate a series from an AR(2) model, holding back the last five
observations to stand in for data that arrives later:

``` r
library(conformalForecast)
#> Loading required package: forecast
set.seed(0)
series_all <- arima.sim(n = 1005, list(ar = c(0.8, -0.5)), sd = 1)
series <- head(series_all, 1000)
new_data <- as.numeric(tail(series_all, 5))
```

Conformal calibration needs a record of past forecast errors. We obtain
one with `cvforecast()`, which refits the model on a rolling forecast
origin and stores the out-of-sample errors at every horizon:

``` r
far2 <- function(x, h, level) {
  Arima(x, order = c(2, 0, 0)) |>
    forecast(h = h, level = level)
}

fc <- cvforecast(series, forecastfun = far2, h = 3, level = c(80, 95),
                 forward = TRUE, initial = 1, window = 100)
```

Now calibrate intervals with split conformal prediction, on a rolling
window of 100 calibration scores:

``` r
scpfc <- conformal(fc, method = "scp", ncal = 100, rolling = TRUE)
```

When the five held-back observations arrive, `update()` extends the fit
incrementally. Intervals that were already issued do not change:

``` r
scpfc_updated <- update(scpfc, new_data = new_data, forecastfun = far2)
```

Finally, evaluate the intervals by forecast horizon:

``` r
# Interval forecast accuracy
accuracy(scpfc_updated, byhorizon = TRUE)
#>        Winkler_95  MSIS_95
#> CV h=1   5.002887 4.830049
#> CV h=2   6.523012 6.297515
#> CV h=3   6.630259 6.397788

# Mean coverage
coverage(scpfc_updated, level = 95)
#>       h=1       h=2       h=3 
#> 0.9503106 0.9476961 0.9400749

# Mean and median interval width
width(scpfc_updated, level = 95, includemedian = TRUE)
#> Mean width:
#>      h=1      h=2      h=3 
#> 4.119936 5.371825 5.406426 
#> 
#> Median width:
#>      h=1      h=2      h=3 
#> 4.074733 5.334842 5.326651
```

## Learn more

- [Introduction to
  conformalForecast](https://xqnwang.github.io/conformalForecast/articles/conformalForecast.html):
  conformal methods, their formulas and tuning arguments, with a
  comparison of coverage and interval width.
- [Updating conformal forecasts with new
  data](https://xqnwang.github.io/conformalForecast/articles/update.html):
  how `update()` extends a fitted pipeline as observations arrive, and
  when it is equivalent to refitting.

## License

This package is free and open source software, licensed under GPL-3.
