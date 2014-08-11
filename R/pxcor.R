#' @rdname corcov
#' @export
pxcor <- function(X, Y, Z) {
  ans <- pxcorArma(X = X, Y = Y, Z = Z) 
  colnames(ans) <- colnames(Y)
  rownames(ans) <- colnames(X)
  return(ans)
}
