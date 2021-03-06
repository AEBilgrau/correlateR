% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RCM.R
\name{ICC}
\alias{ICC}
\alias{Psi2Sigma}
\alias{Sigma2Psi}
\title{RCM miscellaneous functions}
\usage{
ICC(nu, p)

Psi2Sigma(Psi, nu)

Sigma2Psi(Sigma, nu)
}
\arguments{
\item{nu}{A numeric of length one giving the degrees of freedom in the RCM.}

\item{p}{A numeric giving the dimension of the space.}

\item{Psi}{A numeric square positive semi-definite matrix. The underlying 
parameter in the RCM.}

\item{Sigma}{A numeric square positive semi-definite matrix. The expected 
covariance matrix in the RCM.}
}
\value{
\code{ICC}: A numeric giving the ICC.

\code{Psi2Sigma}, \code{Sigma2Psi}:
  The converted matrix the same size as \code{Psi} or \code{Sigma}.
}
\description{
Miscellaneous functions for the random covariance model (RCM).
}
\details{
\code{ICC} compute the ICC in the RCM.
  A simple function for computing the intra-class correlation coefficient 
  (ICC). This function simply computes 1 divided by \code{nu - p}.

\code{Psi2Sigma} and \code{Sigma2Psi} provide conversion between Psi 
  and Sigma. Computes the expected covariance matrix from Psi and nu in the
  random covariance model (RCM) and the other way around.
}
\author{
Anders Ellern Bilgrau <anders.ellern.bilgrau (at) gmail.com>
}
