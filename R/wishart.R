#' Wishart and inverse-Wishart distributions
#'
#' Generate data from Wishart and inverse-Wishart distributions using the 
#' Barlett decomposition.
#'
#' @param n The number of realizations.
#' @param sigma A \code{p} by \code{p} symmetric postive-definite matrix.
#' @param nu An integer giving the degrees of freedom.
#' @return \code{rwishart} returns an \code{p} by \code{p} by \code{n} 
#'   3-dimensional array where each slice is a draw from the Wishart 
#'   distribution.
#' @seealso 
#'   \code{\link{rWishart}} from the \pkg{stats}-package.
#' @examples
#' rwishart(n = 2)
#' rWishart(n = 2, df = 10L, Sigma = diag(3))
#' rwishart(n = 3, diag(4), nu = 4)
#' n <- 1000
#' a <- replicate(n,diag(3)) - rwishart(n)/10L
#' hist(a)
#' @export
rwishart <- function(n, sigma = diag(3), nu = 10L) {
  stopifnot(nrow(sigma) == ncol(sigma))
  p <- ncol(sigma)
  if (nu <= p - 1) {
    warning("For the matrix to be invertible n > p - 1 must hold")
  }
  L <- chol(sigma)
  A <- matrix(0, p, p)
  ans <- replicate(n, {
    A[lower.tri(A)] <- rnorm(p*(p-1)/2)
    diag(A) <- sqrt(rchisq(p, nu - seq_len(p) + 1))
    crossprod(L %*% A)
  })
  return(ans)
}


#' @rdname rwishart
#' @param Psi A square symmetric postive-definite matrix.
#' @return For \code{rinvwishart} each slice is a draw from the inverse Wishart 
#'   distribution.
#' @examples
#' rinvwishart(n = 2)
#' @export
rinvwishart <- function(n, Psi = diag(3), nu = 10L) {
  stopifnot(nrow(Psi) == ncol(Psi))
  stopifnot(nu > ncol(Psi) - 1)
  ans <- replicate(n, solve(rwishart(n = 1, sigma = solve(Psi), nu = nu)[,,1]))
  return(ans)
}
