# Interval estimate accuracy measures

Accuracy measures for interval forecasts.

## Usage

``` r
MSIS(
  lower,
  upper,
  actual,
  train,
  level = 95,
  period,
  d = period == 1,
  D = period > 1,
  na.rm = TRUE,
  ...
)

winkler_score(lower, upper, actual, level = 95, na.rm = TRUE, ...)

interval_measures
```

## Format

An object of class `list` of length 2.

## Arguments

- lower:

  A numeric vector of lower bounds of interval forecasts.

- upper:

  A numeric vector of upper bounds of interval forecasts.

- actual:

  A numeric vector of realised values.

- train:

  A numeric vector of responses used to train the model (for scaled
  scores).

- level:

  The nominal level of the forecast interval (e.g., 95 or 0.95).

- period:

  The seasonal period of the data.

- d:

  Should the response model include a first difference?

- D:

  Should the response model include a seasonal difference?

- na.rm:

  If `TRUE`, remove missing values before calculating the measure.

- ...:

  Additional arguments for each measure.

## Value

For `winkler_score` and `MSIS`, returns a single numeric scalar giving
the average interval score (Winkler or mean scaled interval score).

For the exported object `interval_measures`, returns a **named list of
functions** that can be supplied to higher-level accuracy routines.

## Examples

``` r
set.seed(1)
actual <- rnorm(10)
lower  <- actual - runif(10, 0.5, 1)
upper  <- actual + runif(10, 0.5, 1)
train  <- rnorm(50)

# Winkler score at 95%
winkler_score(lower, upper, actual, level = 95)
#> [1] 1.473882

# Mean scaled interval score (needs training data and period)
MSIS(lower, upper, actual, train, level = 95, period = 1)
#> [1] 1.369701
```
