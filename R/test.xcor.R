#' Test for correlation between multiple variables samples
#' 
#' Test for Pearson correlation between paired samples across many variables.
#' 
#' @param X numeric matrix or vector of data values. 
#' @param Y numeric matrix or vector of data values. 
#'   If matrices, \code{X} and \code{Y} must have the same number of rows.
#'   If vectors, \code{X} and \code{Y} mush have the same length.
#' @param alternative specifies the alternative hypothesis.
#' @param conf.level confidence level for the returned confidence interval.
#' @param null The null to test against. Either a sigle number of or a matrix
#'   of size \code{ncol(X)} times \code{ncol(Y)}.
#' @return 
#'   If matrices are given a list of matrices are returned where the \eqn{ij}th
#'   entry correponds to the results of \code{stats::cor.test(X[,i], Y[,j])}.\cr
#'   If vectors are given a vector of the results is given.
#' @author Anders Ellern Bilgrau <abilgrau (at) math.aau.dk>
#' @details The test and confidence interval is based Fisher's Z transform.
#' @seealso \code{\link[stats]{cor.test}}
#' @examples
#' X <- createData(10, 4)
#' Y <- createData(10, 6)
#' test.xcor(X, Y)
#' test.xcor(X[,1], Y[,1])
#' @export test.xcor
test.xcor <- function(X, Y, 
                      alternative = c("two.sided", "less", "greater"),
                      conf.level = 0.95,
                      null = 0) {
  alternative <- match.arg(alternative)
  if (!is.matrix(X)) {
    dim(X) <- c(length(X), 1)
  }
  if (!is.matrix(Y)) {
    dim(Y) <- c(length(Y), 1)
  }
  stopifnot(nrow(X) == nrow(Y))
  xcor.res <- xcor(X, Y)
  tcor  <- atanh(xcor.res)
  tnull <- atanh(null)
  se <- 1/sqrt(nrow(X) - 3)
  estimate <- tcor - tnull
  
  res <- fisher.test.cor(estimate = estimate, mean = null, se = se, 
                         alternative = alternative, conf.level = conf.level)
  ans <- c(list(xcor = xcor.res), res)
  if (identical(dim(xcor.res), c(1,1))) {
    ans <- unlist(ans)
  }
  return(ans)  
}
