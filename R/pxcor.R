#' @rdname corcov
#' @export
pxcor <- function(X, Y, Z, method = c("Unbiased", "ML")) {
  method <- match.arg(method)
  norm_type <- ifelse(method == "Unbiased", 0L, 1L)
  ans <- pxcovArma(X = X, Y = Y, Z = Z, norm_type = norm_type) 
  colnames(ans) <- colnames(Y)
  rownames(ans) <- colnames(X)
  return(ans)
}
