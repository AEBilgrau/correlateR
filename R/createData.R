#' Create toy data
#' 
#' Creates a named matrix of toy data from an \code{m}-dimensional 
#' zero-mean normal distribution with identity covariance matrix.
#' 
#' @param n an integer giving the number of observations/samples.
#' @param m an integer giving the number of variables/features.
#' @param n.na an integer giving the number of \code{NA}s in the output.
#'   Default is \code{0}.
#' @return A \code{n} times \code{m} numeric matrix of observations.
#' @note The \code{NA} are randomly inserted in the output.
#' @examples
#' createData(10, 4)
#' @export
createData <- function(n, m, n.na = 0) {
  X <- replicate(m, rnorm(n))
  dimnames(X) <- list(paste0("obs", 1:n), paste0("var", 1:m))
  if (n.na) {
    stopifnot(n.na >= 0 && n.na <= n*m)
    X[sample(n*m, n)] <- NA
  }
  return(X)
}
