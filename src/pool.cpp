// Only include aux_functions.h which pulls in RcppArmadillo.h,
// RcppEigen.h, and Rcpp.h

#include "auxiliary_functions.h"

//' Pooled covariance from list of scatter matrices
//' 
//' @param S_list A list of scatter matrices.
//' @param nu A numeric vector giving the number of samples corresponding
//'   to each scatter matrix.
//' @param norm_type A integer of length one equalling \code{0} or \code{1}. 
//'   If \code{0} the bias correction is used. If \code{1} the ML estimate is 
//'   used.
//' @return A numeric matrix giving the pooled variance.
//' @author Anders Ellern Bilgrau
//' @keywords internal
// [[Rcpp::export]]
arma::mat pool(const Rcpp::List & S_list,  // List of scatter matrices
               const Rcpp::NumericVector ns,
               const int norm_type = 0) {
  int k = S_list.size();
  const double denum = sum(ns) - k*(1 - norm_type);
  arma::mat S0 = S_list[0];
  for (int i = 1; i < k; ++i) {
    arma::mat Si = S_list[i];
    S0 += Si;
  }
  return (1.0f/denum) * S0;
}


/*** R
ns <- seq_len(5) + 10
S <- lapply(ns, function(n) correlateR:::scatter(createData(m = 10, n = n)))
correlateR:::pool(S, ns)
*/
