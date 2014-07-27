#' Compute the covariance matrix
#'
#' This function offers the interface to compute the variance-covariance  
#' matrix between the columns of \code{X}. If both \code{X} and \code{Y} is 
#' supplied, the cross-covariances between the columns of the matrices are 
#' computed.
#' 
#' @param X A numeric matrix.
#' @param Y A numeric matrix of compatible dimension with the \code{X}, i.e. 
#'   \code{nrow(X)} equals \code{nrow(Y)}. The default value \code{NULL} is 
#'   equivalent to Y = X howver more efficient.
#' @param method A character of length 1. The unbiased estimate divided with
#'   \code{n-1} and wheras ML uses \code{n}.
#' @return If \code{Y} is not supplied the a square symmetric covariance matrix 
#'   of size \code{ncol(X)} is returned. If \code{Y} is given the 
#'   cross-covariance matrix of size \code{ncol(X)} times \code{ncol(Y)} is 
#'   returned.
#' @details Functions like \link{\code{cov}}. The \code{i}'th and \code{j}'th 
#'   entry of the output matrix  is the covariance between \code{X[i, ]} and 
#'   \code{Y[j, ]}. When \code{Y} is not supplied, then the output is equivalent
#'   to \code{Y} = \code{X}.
#' @author Anders Ellern Bilgrau <abilgrau (at) math.aau.dk>
#' @seealso \link{\code{cov}}
#' @examples
#' X <- replicate(2, rnorm(10))
#' dimnames(X) <- list(paste0("obs", 1:nrow(X)), paste0("var", 1:ncol(X)))
#' Y <- replicate(3, rnorm(10))
#' dimnames(Y) <- list(paste0("obs", 1:nrow(Y)), paste0("var", 1:ncol(Y)))
#' 
#' covariance(X)  # ==  covariance(X, X, "Unbiased")
#' covariance(X, method = "ML")
#' covariance(X, Y, method = "ML")
#' covariance(X, Y, method = "Unbiased")
covariance <- function(X, Y = NULL, method = c("Unbiased", "ML")) {
  method <- match.arg(method)
  norm_type <- ifelse(method == "Unbiased", 0L, 1L)
  
  if (is.null(Y)) {
    ans <- covArma(X = X, norm_type = norm_type) 
    colnames(ans) <- 
    rownames(ans) <- colnames(X)
  } else {
    ans <- crosscovArma(X = X, Y = Y, norm_type = norm_type)
    colnames(ans) <- colnames(Y)
    rownames(ans) <- colnames(X)
  }

  return(ans)
}
