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
  Various functions to compute the marginal (or unconditional) covariance 
  (and cross-covariance) estimates. The functions feature both the maximum
  likelihood and the biased corrected estimates.
*/


// Covariance implementation in Rcpp
// [[Rcpp::export]]
Rcpp::NumericMatrix covRcpp(Rcpp::NumericMatrix & X,
                            const int norm_type) {
  
  const int n = X.nrow();
  const int m = X.ncol();
  const int df = n - 1 + norm_type;
  
  // Centering the matrix!
  X = centerNumericMatrix(X);  // Defined in aux_functions

  // Computing the covariance matrix
  Rcpp::NumericMatrix cov(m, m);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
      cov(i,j) = Rcpp::sum(X(Rcpp::_, i)*X(Rcpp::_, j))/df;
      cov(j,i) = cov(i,j);
    }
  }
  
  return cov;
}


// Cross-covariance implementation in Rcpp
// [[Rcpp::export]]
Rcpp::NumericMatrix xcovRcpp(Rcpp::NumericMatrix & X,
                                 Rcpp::NumericMatrix & Y,
                                 const int norm_type) {
  
  const int n = X.nrow();
  const int m_X = X.ncol();
  const int m_Y = Y.ncol();
  const int df = n - 1 + norm_type;
  
  // Centering the matrices
  X = centerNumericMatrix(X);
  Y = centerNumericMatrix(Y);
  
  // Computing the covariance matrix
  Rcpp::NumericMatrix cov(m_X, m_Y);
  for (int i = 0; i < m_X; ++i) {
    for (int j = 0; j < m_Y; ++j) {
      cov(i,j) = Rcpp::sum(X(Rcpp::_, i)*Y(Rcpp::_, j))/df;
    }
  }
  
  return cov;
}


// covariance "Implementation"" in Armadillio
// [[Rcpp::export]]
arma::mat covArma(const arma::mat & X,
                  const int norm_type) {
  return arma::cov(X, norm_type);
}


// Cross-covariance "Implementation"" in Armadillio
// [[Rcpp::export]]
arma::mat xcovArma(const arma::mat & X,
                       const arma::mat & Y,
                       const int norm_type) {
  return arma::cov(X, Y, norm_type);
}


// covariance in Eigen
// [[Rcpp::export]]
Eigen::MatrixXd covEigen(Eigen::Map<Eigen::MatrixXd> & X,
                         const int norm_type = 0) {

    // Computing degrees of freedom
    // n - 1 is the unbiased estimate whereas n is the MLE
    const int df = X.rows() - 1 + norm_type; // Subtract 1 by default

    // Centering the matrix
    X.rowwise() -= X.colwise().mean();  

    // Return the X^T * X is the scatter matrix
    return X.transpose() * X / df;  // could be .adjoint()
}


// Cross-covariance in Eigen
// [[Rcpp::export]]
Eigen::MatrixXd xcovEigen(Eigen::Map<Eigen::MatrixXd> & X,
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

covArma(X, 0)
covRcpp(X, 0)
covEigen(X, 0)
all.equal(covArma(X, 0), covRcpp(X, 0))
all.equal(covArma(X, 1), covRcpp(X, 1))

all(crosscovRcpp(X, Y, 1) == crosscovArma(X, Y, 1))
print(microbenchmark(crosscovArma(X, Y, 0), 
                     crosscovEigen(X, Y, 0),
                     crosscovRcpp(X, Y, 0), times = 10), order = "median")
*/

