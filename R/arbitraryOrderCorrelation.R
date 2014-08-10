#' Partial correlations
#'
#' Compute arbitrary \code{k}th order partial correlations using recursion. 
#' This implementation is very inefficient.
#' 
##' @param X A numeric matrix with observations in rows and variables
##'   in columns.
#' @param S A numeric correlation matrix
#' @param z A integer vector of indicies of the variabels to condition on. I.e.
#'   length(z) is the order of the partial correlations. If \code{z} has 
#'   length 0 the marginal correlations are returned.
#' @return A numeric matrix of the same size as the correlation matrix
#'   with the partial correlations given the variabels indexed by \code{z}.
#'   \code{NaN} are returned in the rows and columns of the conditioned 
#'   variables.
#' @seealso \code{\link{correlation}} \code{\link{correlation}}
#' @examples
#' X <- createData(6, 20)
#' Y <- createData(6, 20)
#' Z <- createData(6, 20)
#' S <- cor(X)
#' arbitraryOrderCorrelation(S, z = c(2, 3, 5))
#' @export
arbitraryOrderCorrelation <- function(S, z) {
  if (length(z) == 0) {   
    return(S)
  } else if (length(z) == 1) {
    firstOrderPartialCorrelation(S, k = z)
  } else {
    pS <- arbitraryOrderCorrelation(S, z = z[-1])
    firstOrderPartialCorrelation(pS, k = z[1])
  }
}



