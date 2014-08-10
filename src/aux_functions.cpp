// Only include RcppArmadillo.h and RcppEigen.h which pulls in Rcpp.h
#include <RcppArmadillo.h>
#include <RcppEigen.h>

#include "aux_functions.h"

double square(const double x) {
   return x * x;
}

Rcpp::NumericMatrix centerNumericMatrix(Rcpp::NumericMatrix & X) {
  const int m = X.ncol();
  for (int j = 0; j < m; ++j) {
    X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
  }
  return X;
}

// Analogous to lm(y ~ X)$residuals
// [[Rcpp::export]]
arma::colvec residual(const arma::mat & X, const arma::colvec & y) {
    // Construct design matrix (add intercept)
    arma::mat ones = arma::ones<arma::mat>(X.n_rows, 1);
    
    // fit model y = X*coef 
    arma::colvec coef = arma::solve(arma::join_horiz(ones, X), y);
    
    // Compute residuals
    //arma::colvec res = y - arma::join_horiz(ones, X)*coef; // Alternative
    arma::colvec res = (y - coef(0)) - X*coef.subvec(1, coef.n_elem - 1);
    return res;
}