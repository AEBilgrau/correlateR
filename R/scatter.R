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
#' @return The \code{p1} times \code{p1} scatter matrix or 
#'   \code{p1} times \code{p2} cross-scatter matrix.
#' @author Anders Ellern Bilgrau <abilgrau (at) math.aau.dk>
#' @examples
#' n <- 11
#' X <- createData(n, 4)
#' Y <- createData(n, 2)
#' correlateR:::scatter(X)
#' (n-1)*cov(X)
#' 
#' correlateR:::xscatter(X, Y)
#' n*xcov(X, Y, method = "ML")  # == (n - 1)*xcov(X, Y)
#' @export
scatter <- function(X) {
  return(crossprod(center(X)))
}

#' @rdname scatter
#' @export
xscatter <- function(X, Y) {
  return(crossprod(center(X), center(Y)))
}

center <- function(X) {
  return(t(t(X) - colMeans(X)))
}
