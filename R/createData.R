#' Create toy data
#' 
#' Creates a named matrix of toy data from an \code{m}-dimensional 
#' zero-mean normal distribution with identity covariance matrix for 
#' convenience.
#' 
#' @param n an non-negative integer giving the number of observations/samples.
#' @param m an non-negative integer giving the number of variables/features.
#' @param n.na an integer from \code{0} to \code{n*m} giving the number of 
#'   \code{NA}s in the output. Default is \code{0}.
#' @return A \code{n} times \code{m} numeric matrix of observations from
#'   zero-mean gaussian random variables.
#' @note The \code{n.na} \code{NA}s are randomly inserted in the output.
#' @examples
#' createData(10, 4)
#' createData(1, 2)
#' 
#' # Also works for the degenerate cases!
#' createData(3, 0)
#' createData(0, 3)
#' @export
createData <- function(n, m, n.na = 0) {
  stopifnot(n >= 0 && m >= 0)
  X <- matrix(rnorm(n*m), n, m,
              dimnames = list(sprintf("obs%03d", seq_len(n)),
                              sprintf("var%03d", seq_len(m))))
  if (n.na) {
    stopifnot(n.na >= 0 && n.na <= n*m)
    X[sample(n*m, n)] <- NA
  }
  return(X)
}
