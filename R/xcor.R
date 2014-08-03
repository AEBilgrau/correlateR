#' @rdname corcov
#' @export
xcor <- function(X, Y) {
  stopifnot(nrow(X) == nrow(Y))
  ans <- xcorArma(X = X, Y = Y) 
  colnames(ans) <- colnames(X)
  rownames(ans) <- colnames(Y)
  return(ans)
}