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

  A numeric vector of lower bounds of the interval forecasts.

- upper:

  A numeric vector of upper bounds of the interval forecasts.

- actual:

  A numeric vector of realised values.

- train:

  A numeric vector of responses used to train the model. Required by the
  scaled scores.

- level:

  The nominal level of the forecast interval (e.g., 95 or 0.95).
  Defaults to `95`.

- period:

  The seasonal period of the data. Required by the scaled scores.

- d:

  Should the response model include a first difference? Defaults to
  `period == 1`.

- D:

  Should the response model include a seasonal difference? Defaults to
  `period > 1`.

- na.rm:

  If `TRUE`, remove missing values before calculating the measure.
  Defaults to `TRUE`.

- ...:

  Additional arguments for each measure.

## Value

For `winkler_score` and `MSIS`, a single numeric scalar giving the
average interval score (Winkler or mean scaled interval score).

For the exported object `interval_measures`, a named list of functions
that can be supplied to higher-level accuracy routines.

## See also

[`point_measures`](https://xqnwang.github.io/conformalForecast/reference/point_measures.md)
for the point counterparts, and
[`accuracy.cvforecast`](https://xqnwang.github.io/conformalForecast/reference/accuracy.cvforecast.md),
which applies these measures to a cross-validation or a conformal
prediction object.

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
