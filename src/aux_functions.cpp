// Only include RcppArmadillo.h and RcppEigen.h which pulls in Rcpp.h
#include <RcppArmadillo.h>
#include <RcppEigen.h>

#include "aux_functions.h"

double square(const double x) {
   return x * x;
}
