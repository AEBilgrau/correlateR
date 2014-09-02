#' Shrinkage covariance estimation
#'
#' This algorithm seeks to find a covariance estimate that (asymptotically) 
#' minimizes the mean-squared error (MSE) obtained by shrinkage as proposed by
#' Leodit-Wolf (LW). It is effectively a interpolation/mix of the sample ML
#' estimate of the covariance matrix, \eqn{S}, and the most well-conditioned 
#' (and naive) estimate \eqn{F = 1/p tr(S) I}.
#' The algorithm seeks a solution to the problem:
#'   \deqn{\min_{\rho}E[ || \Sigma_O - \Sigma ||^2 ]}{
#'         mininimize E[ || Sigma_O - Sigma ||^2 ] w.r.t. rho}
#'   \deqn{\text{s.t.} \Sigma_O = (1 - \rho)S + \rho F}{
#'         s.t. Sigma_O = (1-rho)*S + rho*F}
#' @param X The data matrix of size \code{n} by \code{p}.
#' @param method The method of estimating the optimal interpolating parameter.
#' @return
#'   A \code{p} by \code{p} numeric matrix with two extra attributes giving
#'   the used mixture (\eqn{\rho}{rho}) and the method.
#' @details
#'   The improved estimate using Rao-Blackwell theorem, abbreviated RBLW, is 
#'   also implemented. 
#'   The intepolated \eqn{\rho}{rho} value used is always 
#'   \eqn{\min(\rho, 1)}{min(rho,1)}.
#' @references
#'   Chen, Y., Member, S. S., Wiesel, A., Eldar, Y. C., Member, S. S., & 
#'   Hero, A. O. (2009). Shrinkage Algorithms for MMSE Covariance Estimation, 
#'   58(734), 1–28. Methodology; Computation. 
#'   Retrieved from http://arxiv.org/abs/0907.4698
#' @examples
#' n <- 3
#' X <- createData(n, 5)
#' cov(X)
#' Scov(X, method = "RBLW")
#' Scov(X, method = "LW")
#' @export
Scov <- function(X, method = c("RBLW", "LW", "OAS")) {
  method <- match.arg(method)

  # Auxiliary functions:
  shrinkage <- function(rho, Shat, Fhat) {  # Interpolation
    return((1 - rho)*Shat + rho*Fhat)
  }
  tr <- function(A) {  # Trace operator
    return(sum(diag(A)))
  }
  frobenious <- function(A) {  # Frobenious norm
    return(sqrt(sum(A^2)))  # == sqrt(tr(tcrossprod(A))))
  }

  n <- nrow(X)    # Number of samples
  p <- ncol(X)    # Number of features
  Shat <- cov(X, method = "ML")   # The ML covariance estimate 1/n, not 1/(n-1)
  Fhat <- diag(tr(Shat)/p, nrow = p) # Well-conditioned est., ie "average variance"

  if (method == "LW") {
    rho <-   # LW estimate of the optimal rho
      sum(sapply(seq_len(n), function(i) frobenious(tcrossprod(X[i,]) - Shat)))/
      (n^2 * tr(Shat^2) - tr(Shat)^2/p)
  } else if (method == "RBLW") {
    rho <- # RBLW estimate of the optimal rho
      ((n - 2)/n * tr(Shat^2) + tr(Shat)^2)/
      ((n + 2) * (tr(Shat^2) - tr(Shat)^2/p))
  } else if (method == "OAS") {
    stop("OAS not implemented yet") 
  } else {
    stop("method ", method, " must be 'RBLW' or 'LW'")
  }
  
  rho_star <- min(rho, 1)
  ans <- shrinkage(rho_star, Shat, Fhat)
  attr(ans, "rho") <- c(rho_used = rho_star, rho = rho)
  attr(ans, "method") <- c(method)
  return(ans)
}
