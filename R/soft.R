#' The soft thresholding function
#' 
#' Element-wise application of the soft thresholding function given by
#' \deqn{soft(a, c) = sgn(a)\cdot(|a| - c)_+}{
#'       soft(a, c) = sgn(a)*max(|a| - c, 0)}
#'
#' @param a A \code{numeric} vector or \code{matrix} to apply the soft 
#'   thresholding function on.
#' @param c A \code{numeric} giving the soft threshold.
#' @return A \code{numeric} vector or \code{matrix} the same size as \code{a}
#'    with thresholded values.
#' @examples
#' a <- rnorm(100)
#' soft(a, c = 0.2)
#'   
#' x <- createData(n = 10, m = 4)
#' soft(x, c = 0.7)
#' @export
soft <- function(a, c) {
  return(sign(a)*pmax(abs(a) - c, 0))
}
