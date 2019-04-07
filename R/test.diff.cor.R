#' Test for difference in correlation
#' 
#' This functions tests the hypothesis of no difference in correlations. It uses
#' the Fisher \eqn{Z} transform (\code{atanh}) to test the null hypothesis
#' of no difference in correlations. See details.
#' 
#' @details 
#'   The function uses the Fisher \eqn{Z} transform (\code{atanh}) of
#'   correlations to test that the hypotheses of no difference in correlation.
#'   The computed \eqn{Z}-score is
#'   \deqn{\frac{Z_1 - Z_2}{\sqrt{1/(n_1 - 3) + 1/(n_2 - 3))}}}{
#'              (Z_1 - Z_2)/ sqrt(1/(n_1 - 3) + 1/(n_2 - 3))}
#'   where \eqn{Z_1} and \eqn{Z_2} are the Fisher transformed correlations.
#'   It performs the test for all correlations in the correlation matrix.
#' 
#' @param X1 A \code{numeric} \code{matrix} of observations.
#' @param X2 A \code{numeric} \code{matrix} of observations.
#' @param cor1 A \code{numeric} \code{matrix} of correlation coefficients in 
#'   the first group.  May be omitted if \code{X1} is provided.
#' @param cor2 A \code{numeric} \code{matrix} of correlation coefficients in 
#'   the second group. May be omitted if \code{X2} is provided.
#' @param n1 \code{integer} of length 1. The number of samples in group 1.
#' @param n2 \code{integer} of length 1. The number of samples in group 2.
#' @param alternative The alternative hypothesis.
#' @param conf.level The confidence level used in the computed confidence 
#'   intervals.
#' @param null A matrix of number giving the difference in correlation under 
#'   the null hypothesis.
#' @return A list of matrices or vector containing:
#'   \item{\code{LCL}}{The lower confidence interval limit.}
#'   \item{\code{UCL}}{The upper confidence interval limit.}
#'   \item{\code{z}}{A numeric matrix of Z-scores for the hypothesis.}
#'   \item{\code{p.val}}{A numeric matrix of the P-values.}
#'   with an attribute giving the alternative hypothesis.
#' @references 
#'   \url{http://core.ecu.edu/psyc/wuenschk/docs30/CompareCorrCoeff.pdf}
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau (at) gmail.com>
#' @seealso 
#'   Similar usage to \code{\link[stats]{cor.test}} (but NOT the same thing).\cr
#'   This is a vectorised version of \code{\link{test.diff.cor.single}}.
#' @examples
#' n1 <- 8
#' n2 <- 10
#' X1 <- createData(n = n1, m = 5)
#' X2 <- createData(n = n2, m = 5)
#' 
#' print(cor1 <- cor(X1))
#' print(cor2 <- cor(X2))
#' 
#' test.diff.cor(X1, X2)
#' 
#' # Directly supplied correlation matrices
#' test.diff.cor(cor1 = cor1, cor2 = cor2, n1 = n1, n2 = n2)
#' 
#' test.diff.cor(X1, X2, alternative = "less")
#' @export test.diff.cor
test.diff.cor <- function(X1, X2,
                          cor1 = cor(X1), 
                          cor2 = cor(X2), 
                          n1 = nrow(X1), 
                          n2 = nrow(X2), 
                          alternative = c("two.sided", "less", "greater"),
                          conf.level = 0.95,
                          null = 0) {
  if (!missing(X1) & !missing(X2)) {
    stopifnot(ncol(X1) == ncol(X2))
  }
  alternative <- match.arg(alternative)
  
  # Fisher z transform the correlations (atanh is the Fisher Z transform)
  Z1 <- atanh(cor1)
  Z2 <- atanh(cor2)
  
  # Compute Z.scores for diff. coexpression
  estimate <- Z1 - Z2
  diag(estimate) <- 0
  se <- sqrt(1/(n1 - 3) + 1/(n2 - 3))

  # Perform the Fisher test
  ans <- fisher.test.cor(estimate = estimate, mean = atanh(null), se = se, 
                         alternative = alternative, conf.level = conf.level)
  attr(ans, "alternative") <- alternative
  return(ans)
}



#' Test for difference in correlation
#' 
#' This functions uses the Fisher Z-transform (atanh) to test the null hypothesis
#' of no difference in correlations between x1 and y1 versus x2 and y2.
#' 
#' @param x1 numeric vector, x-values for the first sample.
#' @param y1 numeric vector, y-values for the first sample.
#' @param x2 numeric vector, x-values for the second sample.
#' @param y2 numeric vector, y-values for the second sample.
#' @return A numeric vector giving correlation for each group, 
#'   size-estimate and standard error, confidence intervals and p-values.
#' @details
#' The \code{alternative} argument specifies the alternative hypothesis given 
#' below.
#' \tabular{rcl}{
#'                    \tab   \tab
#' H0: \code{cor(x1, y1) =  cor(x2, y2)} \cr
#' \code{"two.sided"} \tab =>\tab 
#' H1: \code{cor(x1, y1) != cor(x2, y2)} \cr
#' \code{"greater"}   \tab =>\tab 
#' H1: \code{cor(x1, y1) > cor(x2, y2)} \cr
#' \code{"less"}      \tab =>\tab 
#' H1: \code{cor(x1, y1) < cor(x2, y2)} 
#' }
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau (at) gmail.com>
#' @seealso 
#'   Similar usage to \code{\link[stats]{cor.test}} in \code{stats}, however 
#'   not the same! \cr
#'   See \code{\link{test.diff.cor}} for a vectorised version.
#' @examples
#' x1 <- rnorm(100)
#' y1 <- rnorm(100)
#' x2 <- rnorm(110)
#' y2 <- 4*x2 + 0.5*rnorm(110) + 1
#' 
#' plot(x1, y1, xlim = range(x1, x2), ylim = range(y1, y2))
#' abline(lm(y1 ~ x1))
#' points(x2, y2, col = "red")
#' abline(lm(y2 ~ x2), col = "red")
#' 
#' diff.test <- correlateR:::test.diff.cor.single
#' round(data.frame(
#'   two = diff.test(x1, y1, x2, y2, alternative = "two.sided"),
#'   les = diff.test(x1, y1, x2, y2, alternative = "less"),
#'   gre = diff.test(x1, y1, x2, y2, alternative = "greater")), 2)
#'   
#' round(diff.test(x1, y1, x1, y1, alternative = "two.sided"), 3)
#' round(diff.test(x1, y1, x1, y1, alternative = "less"), 3)
#' round(diff.test(x1, y1, x1, y1, alternative = "less"), 3)
#' @keywords internal
test.diff.cor.single <- 
  function(x1, y1, x2, y2,
           alternative = c("two.sided", "less", "greater"),
           conf.level = 0.95) {
  alternative <- match.arg(alternative)
  stopifnot(length(x1) == length(y1))
  stopifnot(length(x2) == length(y2))
  n1 <- length(x1)
  n2 <- length(x2)
  z1 <- atanh(stats::cor(x1, y1))
  z2 <- atanh(stats::cor(x2, y2))
  var1 <- 1/(n1 - 3)
  var2 <- 1/(n2 - 3)
  se <- sqrt(var1 + var2)
  estimate <- z1 - z2
  z <- estimate/se
  
  if (alternative == "two.sided") {
    q <- -qnorm((1 - conf.level)/2)
    lower <- tanh(estimate - q*se)
    upper <- tanh(estimate + q*se)
    p.val <- 2*pnorm(-abs(z))
  } else if (alternative == "greater") {
    q <- -qnorm((1 - conf.level))
    lower <- tanh(estimate - q*se)
    upper <- 1
    p.val <- 1 - pnorm(z)
  } else if (alternative == "less") {
    q <- -qnorm((1 - conf.level))
    lower <- -1
    upper <- tanh(estimate + q*se)
    p.val <- pnorm(z)
  } else {
    stop("alternative not found")
  }
  ans <- c(r1 = tanh(z1), r2 = tanh(z2),
           estimate = estimate, se = se, z = z,
           lower = lower, upper = upper, p.val = p.val)
  attr(ans, "alternative") <- alternative
  return(ans)         
}
