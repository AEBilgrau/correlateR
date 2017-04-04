// Only include aux_functions.h which pulls in RcppArmadillo.h,
// RcppEigen.h, and Rcpp.h

#include "auxiliary_functions.h"


arma::mat rwishartS(const arma::mat L, const double nu, const int p) {
  // const arma::mat L = arma::chol(sigma);
  arma::mat A(p, p, arma::fill::randn);
  A = trimatl(A); 
  arma::vec chisq(p);
  for (int i = 0; i < p; ++i) {
    chisq[i] = sqrt(Rcpp::as<double>(Rcpp::rchisq(1, nu - i)));
  }
  A.diag() = chisq;
  arma::mat LA = L*A;
  return LA.t()*LA;
}

// [[Rcpp::export]]
arma::cube rwishartArma(const int n, 
                        const arma::mat & sigma, 
                        const double nu) {
  const int p = sigma.n_cols;
  const arma::mat L = arma::chol(sigma, "lower");
  arma::cube ans(p, p, n);
  for (int k = 0; k < n; ++k) {
    ans.slice(k) = rwishartS(L, nu, p);
  }
  return ans;
}

// [[Rcpp::export]]
arma::cube rinvwishartArma(const int n, 
                           const arma::mat & psi, 
                           const double nu) {
  const int p = psi.n_cols ;
  const arma::mat L = arma::chol(inv(psi), "lower");
  arma::cube ans(p, p, n);
  for (int k = 0; k < n; ++k) {
    ans.slice(k) = inv(rwishartS(L, nu, p));
  }
  return ans;
}


/*** R 
nu <- 1e32
p <- 10
correlateR:::rwishartS(diag(p), nu, p)/nu
correlateR:::rwishartArma(n = 2, diag(p), nu)/nu
correlateR:::rinvwishartArma(n = 2, diag(p), nu)/nu

Sigma <- diag(4) 
Sigma[rbind(c(1,2),c(2,1))] <- 0.5
nu <- 1e8

microbenchmark::microbenchmark(
  rWishart(10, nu, Sigma),
  rwishart(10, Sigma, nu),
  rinvwishart(10, Sigma, nu))
*/
