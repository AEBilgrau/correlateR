// Only include aux_functions.h which pulls in RcppArmadillo.h,
// RcppEigen.h, and Rcpp.h

#include "auxiliary_functions.h"


// [[Rcpp::export]]
arma::mat pxcovArma(const arma::mat& X, 
                    const arma::mat& Y, 
                    const arma::mat& Z,
                    const int norm_type = 0L) {
  int m_X = X.n_cols, m_Y = Y.n_cols;
  arma::mat ans(m_X, m_Y);
  
  for (int i = 0; i < m_X; ++i) {
    for (int j = 0; j < m_Y; ++j) {
      arma::colvec rx = residual(Z, X.col(i));
      arma::colvec ry = residual(Z, Y.col(j));
      ans(i,j) = arma::as_scalar(arma::cov(rx, ry, norm_type));
    }
  }
  
  return ans;
}
