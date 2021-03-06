% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fisher.test.cor.R
\name{fisher.test.cor}
\alias{fisher.test.cor}
\title{Fisher's correlation test}
\usage{
fisher.test.cor(estimate, mean, se, alternative, conf.level)
}
\arguments{
\item{estimate}{The fisher z-transformed correlation.}

\item{mean}{The fisher z-transformed null hypotheses.}

\item{se}{The standard error of the estimate (usually \eqn{1/sqrt(n - 3)}
where n is the sample size)}

\item{alternative}{specifies the alternative hypotheses.}

\item{conf.level}{confidence level for the returned confidence interval.}
}
\value{
A list of values resulting from the test.
}
\description{
Test hypotheses about the correlation using fisher's z-transform 
(\code{atanh}).
}
\examples{
est <- matrix(3, 2, 1)
mean <- matrix(0, 2, 1)
sd <- matrix(1, 2, 1)
correlateR:::fisher.test.cor(estimate = est, mean = mean, se = sd,
                             alternative = "two.sided", conf.level = 0.95)
}
\author{
Anders Ellern Bilgrau <anders.ellern.bilgrau (at) gmail.com>
}
\keyword{internal}
