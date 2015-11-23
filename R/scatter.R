#' Compute the scatter matrix
#'
#' This functions is a fast computation of the scatter matrix and cross-scatter
#' matrix given by 
#'   \deqn{\sum_{i = 1}^n(x_i-\mu)(x_i-\mu)^T}{XX^T} and \eqn{XY^T}{XY^T} 
#' respectively.
#' 
#' @param X An \code{n} by \code{p1} numeric matrix with observations in the 
#'   rows and variables in the columns.
#' @param Y An \code{n} by \code{p2} as \code{X}.
#' @param center logical. If \code{TRUE}, the data-matrices are centered first.
#' @return The \code{p1} times \code{p1} scatter matrix or 
#'   \code{p1} times \code{p2} cross-scatter matrix.
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau (at) gmail.com>
#' @examples
#' n <- 11
#' X <- createData(n, 4)
#' Y <- createData(n, 2)
#' scatter(X)
#' (n-1)*cov(X)
#' 
#' xscatter(X, Y)
#' n*xcov(X, Y, method = "ML")  # == (n - 1)*xcov(X, Y)
#' @export
scatter <- function(X, center = TRUE) {
  if (center) {
    return(crossprod(center(X)))
  } else {
    return(crossprod(X))
  }
}

#' @rdname scatter
#' @export
xscatter <- function(X, Y, center = TRUE) {
  if (center) {
    return(crossprod(center(X), center(Y)))
  } else {
    return(crossprod(X, Y))
  }
}

center <- function(X) {
  return(t(t(X) - colMeans(X)))
}
