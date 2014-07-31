#ifndef AUX_FUNCTIONS_H
#define AUX_FUNCTIONS_H

#include <RcppArmadillo.h>
#include <RcppEigen.h>

double square(const double x);

// Subtract the rowmeans of a numeric matrix
Rcpp::NumericMatrix centerNumericMatrix(Rcpp::NumericMatrix & X);

#endif // AUX_FUNCTIONS_H
