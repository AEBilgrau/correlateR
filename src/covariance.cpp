// Only include RcppArmadillo.h and RcppEigen.h which pulls in Rcpp.h
#include <RcppArmadillo.h>
#include <RcppEigen.h>

#include "aux_functions.h"

//// [[Rcpp::depends(RcppArmadillo)]] // Uncomment when sourceCpp()ing
//// [[Rcpp::depends(RcppEigen)]]

//using namespace Rcpp;
//using namespace arma;
//using namespace Eigen;

/*
    Marginal covariance (and cross-covariance) functions
*/

// [[Rcpp::export]]
double cube (const double x) {
  return x * square(x);
}

// [[Rcpp::export]]
arma::mat covArma(const arma::mat & X,
                  const int norm_type) {
  return arma::cov(X, norm_type);
}

// [[Rcpp::export]]
arma::mat crosscovArma(const arma::mat & X,
                       const arma::mat & Y,
                       const int norm_type) {
  return arma::cov(X, Y, norm_type);
}

// [[Rcpp::export]]
Eigen::MatrixXd covEigen(Eigen::Map<Eigen::MatrixXd> & X,
                         const int norm_type = 0) {

    // Computing degrees of freedom
    // n - 1 is the unbiased estimate whereas n is the MLE
    const int df = X.rows() - 1 + norm_type; // Subtract 1 by default

    X.rowwise() -= X.colwise().mean();  // Centering

    // Return the X^T * X is the scatter matrix
    return X.transpose() * X / df;
}

// [[Rcpp::export]]
Eigen::MatrixXd crosscovEigen(Eigen::Map<Eigen::MatrixXd> & X,
                              Eigen::Map<Eigen::MatrixXd> & Y,
                              const int norm_type = 0) {

    // Computing degrees of freedom
    // n - 1 is the unbiased estimate whereas n is the MLE
    const int df = X.rows() - 1 + norm_type; // Subtract 1 by default

    // Centering matrices
    Y.rowwise() -= Y.colwise().mean();
    X.rowwise() -= X.colwise().mean();

    // Return the X^T * X is the scatter matrix
    return X.transpose() * Y / df;
}



/*** R
library("microbenchmark")

X <- replicate(10, rnorm(50))
dimnames(X) <- list(paste0("obs", 1:nrow(X)), paste0("dim", 1:ncol(X)))
Y <- replicate(15, rnorm(50))
dimnames(Y) <- list(paste0("obs", 1:nrow(Y)), paste0("dim", 1:ncol(Y)))

# Covariance check
XX <- unname(stats::cov(X))
stopifnot(all.equal(XX, covArma(X,0)),
          all.equal(XX, crosscovArma(X,X,0)),
          all.equal(XX, covEigen(X,0)),
          all.equal(XX, crosscovEigen(X,X,0)))

microbenchmark(stats::cov(X),
               covArma(X, 0),
               crosscovArma(X, X, 0),
               covEigen(X, 0),
               crosscovEigen(X, X, 0),
               times = 100)


*/

