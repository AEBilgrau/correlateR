#' @rdname pcor
#' @export
pcov <- function(X, z) {
  xpcor(X, X, X[,z])
}