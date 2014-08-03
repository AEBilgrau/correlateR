#' Compute the covariance matrix
#'
#' This function offers the interface to compute the variance-covariance  
#' matrix between the columns of \code{X}. If both \code{X} and \code{Y} is 
#' supplied, the cross-covariances between the columns of the matrices are 
#' computed.
#' 
#' @rdname cor
#' @param method A character of length 1. The unbiased estimate divided with
#'   \code{n-1} and wheras ML uses \code{n}.
#' @return 
#'   \code{cov}: If \code{Y} is not supplied the a square symmetric covariance
#'   matrix of size \code{ncol(X)} is returned. If \code{Y} is given the 
#'   cross-covariance matrix of size \code{ncol(X)} times \code{ncol(Y)} is 
#'   returned.
#' @export
cov <- function(X, method = c("Unbiased", "ML")) {
  method <- match.arg(method)
  norm_type <- ifelse(method == "Unbiased", 0L, 1L)
  ans <- covArma(X = X, norm_type = norm_type) 
  colnames(ans) <- 
  rownames(ans) <- colnames(X)
  return(ans)
}
