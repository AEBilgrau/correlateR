#' @rdname corcov
#' @export
xcov <- function(X, Y, method = c("Unbiased", "ML")) {
  stopifnot(nrow(X) == nrow(Y))
  method    <- match.arg(method)
  norm_type <- ifelse(method == "Unbiased", 0L, 1L)
  ans       <- xcovArma(X = X, Y = Y, norm_type = norm_type) 
  colnames(ans) <- colnames(X)
  rownames(ans) <- colnames(Y)
  return(ans)
}