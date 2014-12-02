#' Wishart and inverse-Wishart distributions
#'
#' Generate data from Wishart and inverse-Wishart distributions using the 
#' Barlett decomposition.
#'
#' @param n The number of realizations.
#' @param sigma A \code{p} by \code{p} symmetric postive-definite matrix.
#' @param nu An numeric of length one giving the degrees of freedom.
#' @return \code{rwishart} returns an \code{p} by \code{p} by \code{n} 
#'   3-dimensional array where each slice is a draw from the Wishart 
#'   distribution.
#' @seealso 
#'   \code{\link{rWishart}} from the \pkg{stats}-package.
#' @examples
#' rwishart(n = 2)
#' rwishart(n = 3, diag(4), nu = 4)
#' @export
rwishart <- function(n, sigma = diag(3), nu = 10) {
  stopifnot(nrow(sigma) == ncol(sigma))
  stopifnot(nu > ncol(sigma) - 1)
  return(rwishartArma(n, sigma, nu))
}

#' @rdname rwishart
#' @param psi A square symmetric postive-definite matrix.
#' @return For \code{rinvwishart} each slice is a draw from the inverse Wishart 
#'   distribution.
#' @examples
#' rinvwishart(n = 2)
#' @export
rinvwishart <- function(n, psi = diag(3), nu = 10) {
  stopifnot(nrow(psi) == ncol(psi))
  stopifnot(nu > ncol(psi) - 1)
  return(rinvwishartArma(n, psi, nu))
}
