#' Partial cross covariance and correlation
#' 
#' Compute the partial cross (x) correlation or covariance using linear 
#' regression
#' 
#' @aliases pxcor pxcov xcov xcor
#' @param x A numeric vector of observations.
#' @param y A numeric vector of observations.
#' @param Z A numeric matrix of observations to condition on.
#' @return \code{pxcor} retuns a numeric  matrix of partial cross correlations
#'   between \code{X} and \code{Y} given \code{Z}.
#' library(RcppEigen)
#' n <- 40
#' X <- replicate(4, rnorm(n))
#' Y <- replicate(2, rnorm(n))
#' Z <- replicate(11, rnorm(n))
#' xcor(X, Y, Z)
#' @export
pxcor <- function(X, Y, Z) {
  pxcorSimple <- function(x, y, Z) {
    Z1 <- cbind(1, Z)
    r_xz <- fastLmPure(Z1, x)$residuals # Alternative: residuals(lm(x ~ Z))
    r_yz <- fastLmPure(Z1, y)$residuals
    return(cor(r_xz, r_yz))
  }
  
  ans <- matrix(NA, ncol(X), ncol(Y))
  for (i in 1:ncol(X)) {
    for (j in 1:ncol(Y)) {
      ans[i,j] <- pxcorSimple(X[, i], Y[, j], Z = Z)
    }
  }
  return(ans)
}


