#' @rdname corcov
#' @param z A integer vector of indices to condition on.
#' @export
pcor <- function(X, z) {
  xpcor(X, X, X[,z])
}