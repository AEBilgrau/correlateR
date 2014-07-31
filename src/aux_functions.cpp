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
