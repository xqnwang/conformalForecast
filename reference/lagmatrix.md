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

  A vector of finite integer lags (positive values) or leads (negative
  values) with a length equal to the number of columns of `x`.

## Value

A matrix with the same class and size as `x`.

## Examples

``` r
set.seed(1)
x <- matrix(rnorm(20), nrow = 5, ncol = 4)

# Create lags of a matrix
lagmatrix(x, c(0, 1, 2, 3))
#>            [,1]       [,2]       [,3]        [,4]
#> [1,] -0.6264538         NA         NA          NA
#> [2,]  0.1836433 -0.8204684         NA          NA
#> [3,] -0.8356286  0.4874291  1.5117812          NA
#> [4,]  1.5952808  0.7383247  0.3898432 -0.04493361
#> [5,]  0.3295078  0.5757814 -0.6212406 -0.01619026
#> attr(,"class")
#> [1] "matrix" "array" 

# Create leads of a matrix
lagmatrix(x, c(0, -1, -2, -3))
#>            [,1]       [,2]       [,3]      [,4]
#> [1,] -0.6264538  0.4874291 -0.6212406 0.8212212
#> [2,]  0.1836433  0.7383247 -2.2146999 0.5939013
#> [3,] -0.8356286  0.5757814  1.1249309        NA
#> [4,]  1.5952808 -0.3053884         NA        NA
#> [5,]  0.3295078         NA         NA        NA
#> attr(,"class")
#> [1] "matrix" "array" 
```
