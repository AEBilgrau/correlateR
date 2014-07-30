// Only include RcppArmadillo.h and RcppEigen.h which pulls in Rcpp.h
#include <RcppArmadillo.h>
#include <RcppEigen.h>

#include "aux_functions.h"

/*
  Marginal correlation (and cross-correlation) functions
*/

// [[Rcpp::export]]
arma::mat corArma(const arma::mat & X,
                  const int norm_type) {
  return arma::cor(X, norm_type);
}

// [[Rcpp::export]]
arma::mat crosscorArma(const arma::mat & X,
                       const arma::mat & Y,
                       const int norm_type) {
  return arma::cor(X, Y, norm_type);
}

// [[Rcpp::export]]
Eigen::MatrixXd corEigen(Eigen::Map<Eigen::MatrixXd> & X,
                         const int norm_type = 0) {
  
  // Computing degrees of freedom
  // n - 1 is the unbiased estimate whereas n is the MLE
  const int df = X.rows() - 1 + norm_type; // Subtract 1 by default
  
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
                              Eigen::Map<Eigen::MatrixXd> & Y,
                              const int norm_type = 0) {
  
  // Computing degrees of freedom
  // n - 1 is the unbiased estimate whereas n is the MLE
  const int df = X.rows() - 1 + norm_type; // Subtract 1 by default
  
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
