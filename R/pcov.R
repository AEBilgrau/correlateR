#' @rdname corcov
#' @export
pcov <- function(X, z) {
  pcovArma(X, z)  # pxcov(X, X, X[, z, drop = FALSE])
}