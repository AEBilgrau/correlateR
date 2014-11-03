#' Fisher's correlation test
#'
#' Test hypothesis about the correlation using fisher's z-transform (atanh).
#' 
#' @param estimate The fisher z-transformed correlation.
#' @param mean The fisher z-transformed null hypothess.
#' @param se The standard error of the estimate (usually \eqn{1/sqrt(n - 3)}
#'   where n is the samplesize)
#' @param alternative specifies the alternative hypothesis.
#' @param conf.level confidence level for the returned confidence interval.
#' @return A list of values resulting from the test.
#' @author Anders Ellern Bilgrau <abilgrau (at) math.aau.dk>
#' @examples
#' est <- matrix(3, 2, 1)
#' mean <- matrix(0, 2, 1)
#' sd <- matrix(1, 2, 1)
#' correlateR:::fisher.test.cor(estimate = est, mean = mean, se = sd,
#'                              alternative = "two.sided", conf.level = 0.95)
#' @keywords internal
fisher.test.cor <- function(estimate, mean, se, alternative, conf.level) {
  z <- (estimate - mean)/se
  l <- drop(matrix(1, nrow(estimate), ncol(estimate)))
  if (alternative == "two.sided") {
    q <- -qnorm((1 - conf.level)/2)
    LCL <- tanh(estimate - q*se)
    UCL <- tanh(estimate + q*se)
    p.val <- 2*pnorm(-abs(z))
  } else if (alternative == "greater") {
    q <- -qnorm((1 - conf.level))
    LCL <- tanh(estimate - q*se)
    UCL <- l
    p.val <- 1 - pnorm(z)
  } else if (alternative == "less") {
    q <- -qnorm((1 - conf.level))
    LCL <- -l
    UCL <- tanh(estimate + q*se)
    p.val <- pnorm(z)
  } else {
    stop("alternative not found")
  }
  return(list(LCL = LCL, UCL = UCL, z = z, p.val = p.val))
}
