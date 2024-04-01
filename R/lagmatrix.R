#' Create lags or leads of a matrix
#'
#' Find a shifted version of a matrix, adjusting the time base backward (lagged)
#' or forward (leading) by a specified number of observations for each column.
#'
#' @param x A matrix or multivariate time series.
#' @param lag A vector of lags (positive values) or leads (negative values) with
#' a length equal to the number of columns of \code{x}.
#'
#' @return A matrix with the same class and size as \code{x}.
#'
#' @examples
#' x <- matrix(rnorm(20), nrow = 5, ncol = 4)
#'
#' # Create lags of a matrix
#' lagmatrix(x, c(0, 1, 2, 3))
#'
#' # Create leads of a matrix
#' lagmatrix(x, c(0, -1, -2, -3))
#'
#' @export
lagmatrix <- function(x, lag) {
  # Ensure 'x' is a matrix
  if (!is.matrix(x))
    stop("ensure x is a matrix")
  n <- nrow(x)
  k <- length(lag)

  if (ncol(x) != k)
    stop("lag must have the same number of columns as x")

  lmat <- x
  for (i in 1:k) {
    if (lag[i] == 0) {
      lmat[, i] <- x[, i]
    } else if (lag[i] > 0) {
      lmat[, i] <- c(rep(NA, lag[i]), x[1:(n - lag[i]), i])
    } else {
      lmat[, i] <- c(x[(abs(lag[i])+1):n, i], rep(NA, abs(lag[i])))
    }
  }
  return(structure(lmat, class = class(x)))
}
