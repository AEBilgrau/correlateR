// Only include aux_functions.h which pulls in RcppArmadillo.h,
// RcppEigen.h, and Rcpp.h

#include "auxiliary_functions.h"

//' Convert covariance matrix to correlation
//' 
//' This functions converts a covariance matrix \code{S} to a correlation matrix
//' fast and efficiently.
//' 
//' @rdname cov2cor
//' @aliases cov2cor
//' @param S A square covariance matrix.
//' @return A square correlation matrix.
//' @author Anders Ellern Bilgrau <abilgrau (at) math.aau.dk>
//' @examples
//' X <- createData(n = 11, m = 4)
//' S <- cov(X)
//' stats::cov2cor(S)
//' cov2corArma(S)
//' if (require(microbenchmark)) {
//'   microbenchmark(A = cov2corArma(S),
//'                  B = stats::cov2cor(S),
//'                  C = cov2cor(S))
//' }
//' @export
// [[Rcpp::export]] 
arma::mat cov2corArma(arma::mat S) {
  arma::colvec inv_sqrt_ss = 1/sqrt(S.diag());
  S.each_col() %= inv_sqrt_ss; // Element-wise multiplication
  S.each_row() %= inv_sqrt_ss.t();
  return S;
}
