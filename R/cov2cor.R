#' @rdname cov2cor
#' @export
cov2cor <- function(S) {
  R <- cov2corArma(S)
  dimnames(R) <- dimnames(S)
  return(R)
}