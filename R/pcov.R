#' @rdname corcov
#' @export
pcov <- function(X, z) {
  pxcor(X, X, X[,z])
}