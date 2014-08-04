#' @rdname corcov
#' @param z A integer vector of indices to condition on.
#' @export
pcor <- function(X, z) {
  pxcor(X, X, X[,z])
}