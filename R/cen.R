#' Centering matrix
#'
#' Constructs a centering matrix of dimension \code{d} by \code{d}.
#'
#' @details Using the centering matrix to center a matrix
#' is not computationally efficient.
#' 
#' @param d integer giving the dimension of the matrix.
#'
#' @return A \code{numeric matrix} of dimension \code{d} by \code{d}.
#' @examples
#' cen(1)
#' cen(2)
#' cen(3L)
#' cen(5)
#' cen(0) # Degenerate case works too
#' @export
cen <- function(d) {
  stopifnot(length(d) == 1)
  if (!isTRUE(all.equal(d, as.integer(d)))) {
    stop("d should be integer-like")
  } 
  return(matrix(-1/d, d, d) + diag(d))
}
