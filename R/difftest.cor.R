#' Test for differential correlation
#' 
#' This functions tests the hypothesis of no difference in correlations.
#' 
#' @param X1 A numeric matrix of observations.
#' @param X2 A numeric matrix of observations.
#' @param cor1 A numeric matrix of correlation coefficients in the first group.
#' @param cor2 As c\code{cor1} for the second group.
#' @param n1 integer of lenth 1. The number of samples in group 1.
#' @param n2 integer of lenth 1. The number of samples in group 2.
#' @param alternative The alternative hypothesis.
#' @param conf.level The confidence level used in the computed confidence 
#'   intervals.
#' @param null A matrix of number giving the difference in correlation under 
#'   the null hypothesis.
#' @return A list of matrices or vector containing:
#'   \item{\code{Z.scores}}{A numeric matrix of Z-scores for the hypothesis.}
#'   \item{\code{P.values}}{A numeric matrix of the P-values.}
#'   \item{\code{lower}}{The lower CI limit.}
#'   \item{\code{upper}}{The upper CI limit.}
#'   with an attribute giving the alternative hypothesis.
#' @details 
#'   The function uses Fisher's Z transform (atanh) of correlations to test
#'   that the hypothes of no difference in correlation. The computed 
#'   Z-score is 
#'   \deqn{\frac{Z1 - Z2}{\sqrt{1/(n1 - 3) + 1/(n2 - 3))}}}{
#'              (Z1 - Z2)/ sqrt(1/(n1 - 3) + 1/(n2 - 3))}
#'   where Z1 and Z2 are the Fisher transformed correlations.
#' @references 
#'   \url{http://core.ecu.edu/psyc/wuenschk/docs30/CompareCorrCoeff.pdf}
#' @author Anders Ellern Bilgrau <abilgrau (at) math.aau.dk>
#' @seealso \code{\link[stats]{cor.test}}
#' @examples
#' n1 <- 8
#' n2 <- 10
#' X1 <- createData(n = n1, m = 5)
#' X2 <- createData(n = n2, m = 5)
#' 
#' print(cor1 <- cor(X1))
#' print(cor2 <- cor(X2))
#' 
#' difftest.cor(X1, X2)
#' difftest.cor(cor1 = cor1, cor2 = cor2, n1 = n1, n2 = n2)
#' 
#' difftest.cor(X1, X2, alternative = "less")
#' @export difftest.cor
difftest.cor <- function(X1, X2,
                         cor1 = cor(X1), cor2 = cor(X2), 
                         n1 = nrow(X1), n2 = nrow(X2), 
                         alternative = c("two.sided", "less", "greater"),
                         conf.level = 0.95,
                         null = 0) {
  if (!missing(X1) & !missing(X2)) {
    stopifnot(ncol(X1) == ncol(X2))
  }
  alternative <- match.arg(alternative)
  
  # Fisher Z transform the correlations (atanh is Fisher's z transformation)
  Z1 <- atanh(cor1)
  Z2 <- atanh(cor2)
  
  # Compute Z.scores for diff. coexpression
  estimate <- Z1 - Z2
  diag(estimate) <- 0
  se <- sqrt(1/(n1 - 3) + 1/(n2 - 3))

  # Perform Fisher's test
  ans <- fisher.test.cor(estimate = estimate, mean = atanh(null), se = se, 
                         alternative = alternative, conf.level = conf.level)
  attr(ans, "alternative") <- alternative
  return(ans)
}

#' Test for difference in correlation
#' 
#' This functions uses Fisher's Z-transform (atanh) to test the null hypothesis
#' of no difference in correlations between x1 and y1 versus x2 and y2.
#' 
#' @param x1 numeric vector, x-values for the first sample.
#' @param y1 numeric vector, y-values for the first sample.
#' @param x2 numeric vector, x-values for the second sample.
#' @param y2 numeric vector, y-values for the second sample.
#' @return A numeric vector giving correlation for each group, 
#'   size-estimate and standard error, confidence intervals and p-values.
#' @details
#'              H0: cor(x1, y1) =  cor(x2, y2) \cr
#' two.sided => H1: cor(x1, y1) != cor(x2, y2) (i.e., z1 - z2 1= 0) \cr
#' greater   => H1: cor(x1, y1) >  cor(x2, y2) (i.e., z1 - z2 > 0) \cr
#' less      => H1: cor(x1, y1) <  cor(x2, y2) (i.e., z1 - z2 < 0)
#' @author Anders Ellern Bilgrau <abilgrau (at) math.aau.dk>
#' @seealso \code{\link[stats]{cor.test}} in \code{stats}
#' @examples
#' x1 <- rnorm(100)
#' y1 <- rnorm(100)
#' x2 <- rnorm(110)
#' y2 <- 4*x2 + 0.5*rnorm(110)
#' 
#' round(data.frame(
#'   two = correlateR:::difftest2.cor(x1, y1, x2, y2),
#'   les =  correlateR:::difftest2.cor(x1, y1, x2, y2, alternative = "less"),
#'   gre = correlateR:::difftest2.cor(x1, y1, x2, y2, alternative = "greater"),
#'   two2 = correlateR:::difftest2.cor(x2, y2, x1, y1)
#'   ),2)
#' @keywords internal
difftest2.cor <- function(x1, y1, x2, y2,
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
