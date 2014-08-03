#' @rdname corcov
#' @export
cov <- function(X, method = c("Unbiased", "ML")) {
  method <- match.arg(method)
  norm_type <- ifelse(method == "Unbiased", 0L, 1L)
  ans <- covArma(X = X, norm_type = norm_type) 
  colnames(ans) <- rownames(ans) <- colnames(X)
  return(ans)
}
