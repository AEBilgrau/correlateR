#' @rdname corcov
#' @export
pcov <- function(X, z) {
  xpcor(X, X, X[,z])
}