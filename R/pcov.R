#' @rdname corcov
#' @export
pcov <- function(X, z) {
  pxcov(X, X, X[, z, drop = FALSE])
}