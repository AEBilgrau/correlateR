#' @rdname corcov
#' @export
pcov <- function(X, z) {
  pcorArma(X, z)  # pxcov(X, X, X[, z, drop = FALSE])
}