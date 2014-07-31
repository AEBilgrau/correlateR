#' Compute the correlation matrix
#'
#' This function offers the interface to compute the correlation  
#' matrix between the columns of \code{X}. If both \code{X} and \code{Y} is 
#' supplied, the cross-correlations between the columns of the matrices are 
#' computed.
#' 
#' @param X A numeric matrix.
#' @param Y A numeric matrix of compatible dimension with the \code{X}, i.e. 
#'   \code{nrow(X)} equals \code{nrow(Y)}. The default value \code{NULL} is 
#'   equivalent to Y = X howver more efficient.
#' @return If \code{Y} is not supplied the a square symmetric correlation matrix 
#'   of size \code{ncol(X)} is returned. If \code{Y} is given the 
#'   cross-correlation matrix of size \code{ncol(X)} times \code{ncol(Y)} is 
#'   returned.
#' @details Functions like \code{\link{cov}}. The \code{i}'th and \code{j}'th 
#'   entry of the output matrix  is the correlation between \code{X[i, ]} and 
#'   \code{Y[j, ]}. When \code{Y} is not supplied, then the output is equivalent
#'   to \code{Y} = \code{X}.
#' @author Anders Ellern Bilgrau <abilgrau (at) math.aau.dk>
#' @seealso \code{\link{cov}}
#' @examples
#' X <- replicate(2, rnorm(10))
#' dimnames(X) <- list(paste0("obs", 1:nrow(X)), paste0("var", 1:ncol(X)))
#' Y <- replicate(3, rnorm(10))
#' dimnames(Y) <- list(paste0("obs", 1:nrow(Y)), paste0("var", 1:ncol(Y)))
#' 
#' correlation(X)
#' correlation(X, X)
#' correlation(X, Y)
correlation <- function(X, Y = NULL) {
  if (is.null(Y)) {
    ans <- corArma(X = X) 
    colnames(ans) <- 
    rownames(ans) <- colnames(X)
  } else {
    stopifnot(nrow(X) == nrow(Y))
    ans <- crosscorArma(X = X, Y = Y)
    colnames(ans) <- colnames(Y)
    rownames(ans) <- colnames(X)
  }
  return(ans)
}
