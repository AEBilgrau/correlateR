#' Shrinkage covariance estimation
#'
#' This algorithm seeks to find a covariance (dense) estimate
#' that (asymptotically) minimizes the mean-squared error (MSE) obtained by
#' linear shrinkage problem as proposed by Ledoit and Wolf (LW).
#' It is effectively a interpolation/mix of the sample ML
#' estimate of the covariance matrix, \eqn{S}, and the most well-conditioned
#' (and naive) estimate \eqn{F = 1/p tr(S) I}.
#'
#' @details
#'   The improved estimate using Rao-Blackwell theorem, abbreviated RBLW, and
#'   the oracle approximating shrinkage (OAS) are also implemented. The 
#'   algorithm seeks a solution to the problem:
#'     \deqn{\min_{\rho}E[ || \Sigma_O - \Sigma ||^2 ]}{%
#'           mininimize E[ || Sigma_O - Sigma ||^2 ] w.r.t. rho}
#'     \deqn{s.t. \Sigma_O = (1 - \rho)S + \rho F}{%
#'           s.t. Sigma_O = (1-rho)*S + rho*F}
#'   using various methods
#'   The intepolated \eqn{\rho}{rho} value used is always
#'   \eqn{\min(\rho, 1)}{min(rho,1)}.
#'   More information can be found in the given reference.
#'
#' @param X The data matrix of size \code{n} by \code{p}.
#' @param method The method of estimating the optimal interpolating parameter.
#'   The default is OAS.
#' @return
#'   A \code{p} by \code{p} numeric matrix with two extra attributes giving
#'   the used mixture (\eqn{\rho}{rho}) and the method.
#' @references
#'   Ledoit, O., & Wolf, M. (2004). A well-conditioned estimator for
#'   large-dimensional covariance matrices. Journal of Multivariate Analysis,
#'   88(2), 365-411. doi:10.1016/S0047-259X(03)00096-4
#'
#'   Chen, Y., & Wiesel, A. (2010). Shrinkage algorithms for MMSE covariance
#'   estimation. Signal Processing, IEEE, 58(734), 1-28.
#'   Methodology; Computation.
#'   \url{http://arxiv.org/abs/0907.4698}
#'
#'   Schafer, J., & Strimmer, K. (2005). A shrinkage approach to large-scale
#'   covariance matrix estimation and implications for functional genomics.
#'   Statistical Applications in Genetics and Molecular Biology, 4(1).
#'   \url{http://www.stat.wisc.edu/courses/st992-newton/smmb/files/expression/shrinkcov2005.pdf}
#' @examples
#' n <- 3
#' X <- createData(n, 5)
#' cov(X)
#' Scov(X, method = "OAS")
#' Scov(X, method = "RBLW")
#' Scov(X, method = "LW")
#' @aliases Scov Scor
#' @export
Scov <- function(X, method = c("OAS", "RBLW", "LW", "SS")) {
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
    rho <-   # LW estimate of the optimal rho (eq. 13 i ref.)
      sum(sapply(seq_len(n), function(i) frobenious(tcrossprod(X[i,]) - Shat)))/
      (n^2 * tr(Shat^2) - tr(Shat)^2/p)
  } else if (method == "RBLW") {
    rho <- # RBLW estimate of the optimal rho (eq. 17 i ref.)
      ((n - 2)/n * tr(Shat^2) + tr(Shat)^2)/
      ((n + 2) * (tr(Shat^2) - tr(Shat)^2/p))
  } else if (method == "OAS") {
    rho <- # OAS estimate of the optimal rho (eq. 23 i ref.)
      ((1 - 2/p) * tr(Shat^2) + tr(Shat)^2)/
      ((n + 1 - 2/p) * (tr(Shat^2) - tr(Shat)^2/p))
  } else {
    stop("method ", method, " must one 'OAS', 'RBLW', 'LW'. ",
         "Others (SS) are not implemented yet.")
  }

  rho_star <- min(rho, 1)
  ans <- shrinkage(rho_star, Shat, Fhat)
  attr(ans, "rho") <- c("used.rho" = rho_star, "rho" = rho)
  attr(ans, "method") <- c(method)
  return(ans)
}
