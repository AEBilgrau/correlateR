#' Wishart and inverse-Wishart distributions
#'
#' Generate data from a Wishart and inverse-Wishart distribution.
#'
#' @param n The number of realizations.
#' @param sigma A square symmetric postive-definite matrix.
#' @param nu An integer giving the degrees of freedom.
#' @return \code{rwishart} returns a matrix drawn from the Wishart distribution.
#' @seealso 
#'   \code{\link{rWishart}} from the \pkg{stats}-package.
#' @examples
#' rwishart(n = 2)
#' rWishart(n = 2, df = 10L, Sigma = diag(3))
#' rwishart(n = 3, diag(4), nu = 4)
#' n <- 10000
#' a <- replicate(n,diag(3)) - rwishart(n)/10L
#' hist(a)
#' @export
rwishart <- function(n, sigma = diag(3), nu = 10L) {
  stopifnot(nrow(sigma) == ncol(sigma))
  p <- ncol(sigma)
  if (nu <= p - 1) {
    warning("For the matrix to be invertible n > p - 1 must hold")
  }
  ans <- replicate(n, crossprod(rmvnormal(nu, mu = rep(0, p), sigma = sigma)))
  return(ans)
}

#' @rdname rwishart
#' @param Psi A square symmetric postive-definite matrix.
#' @return \code{rinvwishart} returns a matrix drawn from the inverse
#'   Wishart distribution.
#' @examples
#' rinvwishart(n = 2)
#' @export
rinvwishart <- function(n, Psi = diag(3), nu = 10L) {
  stopifnot(nrow(Psi) == ncol(Psi))
  stopifnot(nu > ncol(Psi) - 1)
  ans <- replicate(n, solve(rwishart(n = 1, sigma = solve(Psi), nu = nu)[,,1]))
  return(ans)
}
