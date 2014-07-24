// Only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]] // Uncomment when sourceCpp()ing
// Use RcppEigen??

//using namespace Rcpp;
//using namespace RcppArmadillo;

// [[Rcpp::export]]
arma::mat covRcpp(arma::mat X) {
  return arma::cov(X);
}

// [[Rcpp::export]]
arma::mat corRcpp(arma::mat X) {
  return arma::cor(X);
}


/*** R
library("microbenchmark")

X <- replicate(1000, rnorm(1000));
dimnames(X) <- list(paste0("obs", 1:nrow(X)), paste0("dim", 1:ncol(X)))

# Covariance check
all.equal(unname(cov(X)), covRcpp(X))
microbenchmark(cov(X), covRcpp(X))

# Correlation check
all.equal(unname(stats::cor(X)), corRcpp(X))
microbenchmark(stats::cov(X), 
               covRcpp(X), 
               WGCNA::corFast(X, nThreads = 1),
               times = 1)

*/

