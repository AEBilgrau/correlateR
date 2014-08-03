#' Compute the correlation matrix
#'
#' This function offers the interface to compute the correlation  
#' matrix between the columns of \code{X}.
#' 
#' @aliases pcor pcov cor cov
#' @param X A numeric matrix with observations in rows and variables/features
#'   in columns.
#' @return 
#'   \code{cor}: If \code{Y} is not supplied the a square symmetric correlation  
#'   matrix of size \code{ncol(X)} is returned. If \code{Y} is given the 
#'   cross-correlation matrix of size \code{ncol(X)} times \code{ncol(Y)} is 
#'   returned.
#' @details 
#'   Functions like \code{\link{cov}}. The \code{i}'th and \code{j}'th 
#'   entry of the output matrix  is the correlation between \code{X[i, ]} and 
#'   \code{Y[j, ]}. When \code{Y} is not supplied, then the output is equivalent
#'   to \code{Y} = \code{X}.
#' @author Anders Ellern Bilgrau <abilgrau (at) math.aau.dk>
#' @seealso \code{\link{corFamily}} (the workhorse functions)
#' @examples
#' X <- replicate(2, rnorm(10))
#' dimnames(X) <- list(paste0("obs", 1:nrow(X)), paste0("var", 1:ncol(X)))
#' 
#' cor(X)
#' cov(X, method = "ML")
#' cov(X, method = "Unbiased")
#' @export
cor <- function(X) {

    ans <- corArma(X = X) 
    colnames(ans) <- 
    rownames(ans) <- colnames(X)

  return(ans)
}
