#' Covariance and correlation
#'
#' These functions offers a interface to compute arbitrary order partial
#' (or semi-partial) variance-covariance or correlation matrix as well as
#' partial cross variance-covariances or correlations.
#'
#' @rdname corcov
#' @aliases cor cov pcor pcov xcor xcov pxcor pxcov correlation covariance
#' @param X A numeric matrix with observations in rows and variables/features
#'   in columns.
#' @param Y A numeric matrix with the same number of rows as \code{X}.
#' @param Z A numeric matrix with the same number of rows as \code{X}. This is
#'   the 
#' @param method A character of length 1. The unbiased estimate divided with
#'   \code{n-1} and wheras ML uses \code{n}.
#' @return 
#'   All functions return a matrix of correlations or covariances.
#' @details 
#'   Give some details!
#' @author Anders Ellern Bilgrau <abilgrau (at) math.aau.dk>
#' @seealso \code{\link{corFamily}} (the workhorse functions)
#' @examples
#' n <- 11
#' X <- createData(n, 4)
#' Y <- createData(n, 2)
#' Z <- createData(n, 9)
#' 
#' cor(X)
#' cov(X, method = "ML")
#' cov(X, method = "Unbiased")
#' 
#' xcov(X, Y)
#' xcor(X, Y)
#' 
#' pcov(X, z = numeric(0))  # == cov(X)
#' pcov(X, z = 1:3)
#' pcor(X, z = numeric(0))
#' pcor(X, z = 3)
#' 
#' pxcov(X, Y, Z)
#' pxcor(X, Y, Z)
#' @export
cor <- function(X) {
  ans <- corArma(X = X) 
  colnames(ans) <- rownames(ans) <- colnames(X)
  return(ans)
}
