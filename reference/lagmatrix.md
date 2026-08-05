# Create lags or leads of a matrix

Find a shifted version of a matrix, adjusting the time base backward
(lagged) or forward (leading) by a specified number of observations for
each column.

## Usage

``` r
lagmatrix(x, lag)
```

## Arguments

- x:

  A matrix or multivariate time series.

- lag:

  A vector of lags (positive values) or leads (negative values) with a
  length equal to the number of columns of `x`.

## Value

A matrix with the same class and size as `x`.

## Examples

``` r
x <- matrix(rnorm(20), nrow = 5, ncol = 4)

# Create lags of a matrix
lagmatrix(x, c(0, 1, 2, 3))
#>            [,1]         [,2]       [,3]      [,4]
#> [1,]  0.4755095           NA         NA        NA
#> [2,] -0.7099464  0.291446236         NA        NA
#> [3,]  0.6107264 -0.443291873 -0.5686687        NA
#> [4,] -0.9340976  0.001105352 -0.1351786 0.3329504
#> [5,] -1.2536334  0.074341324  1.1780870 1.0630998
#> attr(,"class")
#> [1] "matrix" "array" 

# Create leads of a matrix
lagmatrix(x, c(0, -1, -2, -3))
#>            [,1]         [,2]       [,3]      [,4]
#> [1,]  0.4755095 -0.443291873  1.1780870 0.3700188
#> [2,] -0.7099464  0.001105352 -1.5235668 0.2670988
#> [3,]  0.6107264  0.074341324  0.5939462        NA
#> [4,] -0.9340976 -0.589520946         NA        NA
#> [5,] -1.2536334           NA         NA        NA
#> attr(,"class")
#> [1] "matrix" "array" 
```
