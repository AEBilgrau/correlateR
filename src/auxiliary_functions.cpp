// Only include RcppArmadillo.h and RcppEigen.h which pulls in Rcpp.h
#include <RcppArmadillo.h>
#include <RcppEigen.h>

#include "auxiliary_functions.h"

double square(const double x) {
   return x * x;
}

// Function for centering matrix
Rcpp::NumericMatrix centerNumericMatrix(Rcpp::NumericMatrix & X) {
  const int m = X.ncol();
  for (int j = 0; j < m; ++j) {
    X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
  }
  return X;
}

// Analogous to lm(y ~ X)$residuals.
arma::colvec residual(const arma::mat & X, const arma::colvec & y) {
    // Return y if X is degenerate
    if (X.n_rows == 0 || X.n_cols == 0) {
      return y;
    }
    
    // Construct design matrix (add intercept)
    arma::mat ones = arma::ones<arma::mat>(X.n_rows, 1);
    
    // fit model y = X*coef 
    arma::colvec coef = arma::solve(arma::join_horiz(ones, X), y);
    
    // Compute residuals
    //arma::colvec res = y - arma::join_horiz(ones, X)*coef; // Alternative
    arma::colvec res = (y - coef(0)) - X*coef.subvec(1, coef.n_elem - 1);
    return res;
}

// Multivariate (p-dimensional) log gamma function
// See Bmisc for further details
// [[Rcpp::export]]
Rcpp::NumericVector lgammap(const Rcpp::NumericVector & x, const int p = 1) {
  const double c0 = log(M_PI)*p*(p - 1)/4;
  Rcpp::NumericVector ans(x.size(), c0);
  for (int j = 0; j < p; ++j) {
    ans += Rcpp::lgamma(x - j/2.0f);
  }
  return ans;
}


// Compute log determinant
// 
// Compute the sign and log of the determinant. Faster than \code{determinant}.
// 
// @param x A numeric matrix.
// @return Returns a 1 by 2 matrix where the first element is the log of the 
//  determinant and the second element is the sign of the determinant.
// @seealso \code{\link{determinant}}
// @author Anders Ellern Bilgrau
// @examples
// y <- matrix(rnorm(100), 10, 10)
// determinant(y)
// correlateR:::logdet_arma(y)
// [[Rcpp::export]]
arma::vec logdet_arma(const arma::mat & x) {
  arma::vec val_sign(2); // First element value and second element is sign
  log_det(val_sign(0), val_sign(1), x);
  return val_sign;
}


