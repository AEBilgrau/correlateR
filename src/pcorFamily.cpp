// Only include aux_functions.h which pulls in RcppArmadillo.h,
// RcppEigen.h, and Rcpp.h

#include "auxiliary_functions.h"

// [[Rcpp::export]]
arma::mat pcorArma(const arma::mat& X, const arma::uvec& z) {
  
  int m = X.n_cols;
  arma::uvec zz = z - 1;  // Remember, shift index by one!
  arma::mat ans(m, m);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {

      arma::colvec rx = residual(X.cols(zz), X.col(i));
      arma::colvec ry = residual(X.cols(zz), X.col(j));
      ans(i, j) = arma::as_scalar(arma::cor(rx, ry));
      ans(j, i) = ans(i,j);  // Autocorrelation is symmetric
      
    }
  }

  return ans;
}
