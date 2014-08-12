#ifndef AUXILIARY_FUNCTIONS_H
#define AUXILIARY_FUNCTIONS_H

#include <RcppArmadillo.h>
#include <RcppEigen.h>

double square(const double x);

// Subtract the rowmeans of a numeric matrix
Rcpp::NumericMatrix centerNumericMatrix(Rcpp::NumericMatrix & X);

// Compute residuals in linear model
arma::colvec residual(const arma::mat & X, const arma::colvec & y);

#endif // AUXILIARY_FUNCTIONS_H
