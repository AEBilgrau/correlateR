// Only include RcppArmadillo.h and RcppEigen.h which pulls in Rcpp.h
#include <RcppArmadillo.h>
#include <RcppEigen.h>

#include "aux_functions.h"

/*
  Various functions to compute the marginal (or unconditional) correlations 
  (and cross-covariance) estimates. The functions feature both the maximum
  likelihood and the biased corrected estimates.
*/

// Correlation implementation in Rcpp
// [[Rcpp::export]]
Rcpp::NumericMatrix corRcpp(Rcpp::NumericMatrix & X) {
  
  const int m = X.ncol();
  
  // Centering the matrix
  X = centerNumericMatrix(X);
  
  // Compute 1 over the sample standard deviation
  Rcpp::NumericVector inv_sqrt_ss(m);
  for (int i = 0; i < m; ++i) {
    inv_sqrt_ss(i) = 1/sqrt(Rcpp::sum(X(Rcpp::_, i)*X(Rcpp::_, i)));
  }
  
  // Computing the correlation matrix
  Rcpp::NumericMatrix cor(m, m);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < i; ++j) {
      cor(i, j) = Rcpp::sum(X(Rcpp::_,i)*X(Rcpp::_,j)) *
        inv_sqrt_ss(i) * inv_sqrt_ss(j);
      cor(j, i) = cor(i, j);
    }
    cor(i,i) = 1;
  }

  return cor;
}


// Cross-correlation implementation in Rcpp
// [[Rcpp::export]]
Rcpp::NumericMatrix crosscorRcpp(Rcpp::NumericMatrix & X,
                                 Rcpp::NumericMatrix & Y) {
  
  const int m_X = X.ncol();
  const int m_Y = Y.ncol();
  
  // Centering the matrices
  X = centerNumericMatrix(X);
  Y = centerNumericMatrix(Y);
  
  
  // Compute 1 over square root the sum of squares
  Rcpp::NumericVector inv_sqrt_ss_X(m_X);
  for (int i = 0; i < m_X; ++i) {
    inv_sqrt_ss_X(i) = 1/sqrt(Rcpp::sum(X(Rcpp::_, i)*X(Rcpp::_, i)));
  }
  Rcpp::NumericVector inv_sqrt_ss_Y(m_Y);
  for (int i = 0; i < m_Y; ++i) {
    inv_sqrt_ss_Y(i) = 1/sqrt(Rcpp::sum(Y(Rcpp::_, i)*Y(Rcpp::_, i)));
  }
  
  // Computing the cross-correlation matrix
  Rcpp::NumericMatrix cor(m_X, m_Y);
  for (int i = 0; i < m_X; ++i) {
    for (int j = 0; j < m_Y; ++j) {
      cor(i, j) = Rcpp::sum(X(Rcpp::_, i)*Y(Rcpp::_, j)) *
        inv_sqrt_ss_X(i) * inv_sqrt_ss_Y(j);
    }
  }

  return cor;
}



// Cross-correlation implementation in Rcpp
// [[Rcpp::export]]
arma::mat corArma(const arma::mat & X) {
  return arma::cor(X, 0);
}

// Cross-correlation implementation in Rcpp
// [[Rcpp::export]]
arma::mat crosscorArma(const arma::mat & X,
                       const arma::mat & Y) {
  return arma::cor(X, Y, 0);
}



// [[Rcpp::export]]
Eigen::MatrixXd corEigen(Eigen::Map<Eigen::MatrixXd> & X) {
  
  // Computing degrees of freedom
  // n - 1 is the unbiased estimate whereas n is the MLE
  const int df = X.rows() - 1; // Subtract 1 by default
  
  X.rowwise() -= X.colwise().mean();  // Centering
  
  Eigen::MatrixXd cor = X.transpose() * X / df;   // The covariance matrix
  
  // Get 1 over the standard deviations
  Eigen::VectorXd inv_sds = cor.diagonal().array().sqrt().inverse();
  
  // Scale the covariance matrix
  cor = cor.cwiseProduct(inv_sds * inv_sds.transpose());
  
  return cor;
}

// [[Rcpp::export]]
Eigen::MatrixXd crosscorEigen(Eigen::Map<Eigen::MatrixXd> & X,
                              Eigen::Map<Eigen::MatrixXd> & Y) {
  
  // Computing degrees of freedom
  // n - 1 is the unbiased estimate whereas n is the MLE
  const int df = X.rows() - 1; // Subtract 1 by default
  
  // Centering matrices
  Y.rowwise() -= Y.colwise().mean();
  X.rowwise() -= X.colwise().mean();
  
  // The covariance matrix
  Eigen::MatrixXd cor = X.transpose() * Y / df;
  
  // Compute 1 over the standard deviations of X and Y
  Eigen::VectorXd inv_sds_X = (X.colwise().norm()/sqrt(df)).array().inverse();
  Eigen::VectorXd inv_sds_Y = (Y.colwise().norm()/sqrt(df)).array().inverse();
  
  // Scale the covariance matrix
  cor = cor.cwiseProduct(inv_sds_X * inv_sds_Y.transpose());
  return cor;
}


