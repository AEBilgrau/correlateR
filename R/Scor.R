#' @rdname Scov
#' @export
Scor <- function(X, method = c("OAS", "RBLW", "LW", "SS")) {
  method <- match.arg(method)
  return(cov2cor(Scov(X, method = method)))
}
