# Point estimate accuracy measures

Accuracy measures for point forecast residuals.

## Usage

``` r
ME(resid, na.rm = TRUE)

MAE(resid, na.rm = TRUE, ...)

MSE(resid, na.rm = TRUE, ...)

RMSE(resid, na.rm = TRUE, ...)

MPE(resid, actual, na.rm = TRUE, ...)

MAPE(resid, actual, na.rm = TRUE, ...)

MASE(
  resid,
  train,
  demean = FALSE,
  na.rm = TRUE,
  period,
  d = period == 1,
  D = period > 1,
  ...
)

RMSSE(
  resid,
  train,
  demean = FALSE,
  na.rm = TRUE,
  period,
  d = period == 1,
  D = period > 1,
  ...
)

point_measures
```

## Format

An object of class `list` of length 8.

## Arguments

- resid:

  A numeric vector of residuals from either the validation or test data.

- na.rm:

  If `TRUE`, remove missing values before calculating the measure.

- ...:

  Additional arguments for each measure.

- actual:

  A numeric vector of responses matching the forecasts (for percentage
  measures).

- train:

  A numeric vector of responses used to train the model (for scaled
  measures).

- demean:

  Should the response be demeaned (for MASE and RMSSE)?

- period:

  The seasonal period of the data.

- d:

  Should the response model include a first difference?

- D:

  Should the response model include a seasonal difference?

## Value

For the individual functions (`ME`, `MAE`, `MSE`, `RMSE`, `MPE`, `MAPE`,
`MASE`, `RMSSE`), returns a single numeric scalar giving the requested
accuracy measure.

For the exported object `point_measures`, returns a **named list of
functions** that can be supplied to higher-level accuracy routines.

## Examples

``` r
# Toy residuals and data
set.seed(1)
y_train <- rnorm(50)
y_test  <- rnorm(10)
fcast   <- y_test + rnorm(10, sd = 0.2)
resid   <- y_test - fcast

# Basic measures
ME(resid)
#> [1] -0.09024199
MAE(resid)
#> [1] 0.1937409
RMSE(resid)
#> [1] 0.2606406

# Percentage measures require 'actual'
MPE(resid, actual = y_test)
#> [1] 9.271545
MAPE(resid, actual = y_test)
#> [1] 62.86441

# Scaled measures require training data (and seasonal period if applicable)
MASE(resid, train = y_train, period = 1)
#> [1] 0.2124482
RMSSE(resid, train = y_train, period = 1)
#> [1] 0.2271102
```
