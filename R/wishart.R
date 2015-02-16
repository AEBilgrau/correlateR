#' Wishart and inverse-Wishart distributions
#'
#' Generate data from Wishart and inverse-Wishart distributions using the 
#' Barlett decomposition.
#' When simulating from the Wishart distribution and \code{nu} is an integer
#' strictly less than p, then the functions draws from a singular Wishart 
#' distribution using multivariate normal draws.
#'
#' @param n The number of realizations.
#' @param sigma A \code{p} by \code{p} symmetric positive definite matrix.
#' @param nu An numeric of length one giving the degrees of freedom.
#' @param verbose logical. Warn if draws are from a singular Wishart 
#'   distribution.
#' @return \code{rwishart} returns an \code{p} by \code{p} by \code{n} 
#'   3-dimensional array where each slice is a draw from the Wishart 
#'   distribution.
#' @seealso 
#'   \code{\link{rWishart}} from the \pkg{stats}-package.
#' @examples
#' rwishart(n = 2)
#' rwishart(n = 3, diag(4), nu = 4)
#' 
#' rwishart(n = 1, diag(3), nu = 3)
#' rwishart(n = 1, diag(3), nu = 2)
#' rwishart(n = 1, diag(3), nu = 1)
#' @importFrom GMCM rmvnormal
#' @export
rwishart <- function(n, sigma = diag(3), nu = 10, verbose = TRUE) {
  stopifnot(nrow(sigma) == ncol(sigma))
  if (nu <= ncol(sigma) - 1) {
    if (!all.equal(nu, as.integer(nu))) {
      stop("When nu <= p - 1 then nu must be an integer.")
    }
    if (verbose) {
      warning("Simulating from a singular Wishart distribution as nu <= p - 1.")
    }
    mu <- rep(0, nrow(sigma))
    return(replicate(n, crossprod(rmvnormal(nu, mu = mu, sigma = sigma))))
  }
  return(rwishartArma(n, sigma, nu))
}

#' @rdname rwishart
#' @param psi A square symmetric positive definite matrix.
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
