#include <RcppArmadillo.h>
#include <RcppEigen.h>

#include "aux_functions.h"


// [[Rcpp::export]]
arma::mat pcorArma(const arma::mat& X, const arma::uvec& z) {
  int m = X.n_cols;
  arma::mat ans(m, m);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
      // Remember, shift index by one!
      arma::colvec rx = residual(X.cols(z-1u), X.col(i));
      arma::colvec ry = residual(X.cols(z-1u), X.col(j));
      ans(i,j) = arma::as_scalar(arma::cor(rx, ry));
      ans(j,i) = ans(i,j);
    }
  }

  return ans;
}
