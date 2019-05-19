// Only include aux_functions.h which pulls in RcppArmadillo.h,
// RcppEigen.h, and Rcpp.h

#include "auxiliary_functions.h"

//// [[Rcpp::depends(RcppArmadillo)]] // Uncomment when sourceCpp()ing
//// [[Rcpp::depends(RcppEigen)]]

//' Marginal covariance matrix
//' 
//' Various workhorse functions to compute the marginal (or unconditional)
//' covariance (and cross-covariance) estimates. The functions feature both the
//' maximum likelihood and the biased corrected estimates. They are (almost)
//' equivalent implementations of \code{\link[stats:cor]{cov}} (\code{stats::cov})
//' in Rcpp, RcppArmadillo, and RcppEigen.
//' 
//' @rdname covFamily
//' @aliases covFamily covRcpp xcovRcpp covArma xcovArma covEigen xcovEigen
//' @param X A numeric matrix.
//' @param Y A numeric matrix of compatible dimension with the \code{X}, i.e. 
//'   \code{nrow(X)} equals \code{nrow(Y)}.
//' @param norm_type an integer of length one giving the estimator. The 
//'   default \code{0L} gives the unbiased estimate while \code{1L} gives the 
//'   MLE.
//' @return
//'   The \code{cor}-family returns a numeric correlation matrix of size 
//'   \code{ncol(X)} times \code{ncol(X)}.
//'   
//'   The \code{xcor}-family returns a numeric cross-covariance matrix 
//'   of size \code{ncol(X)} times \code{ncol(Y)}.
//' @details
//'   Functions almost like \code{\link{cor}}.
//'   For the \code{xcorXX} functions, the \code{i}'th and \code{j}'th 
//'   entry of the output matrix is the correlation between \code{X[i, ]} and 
//'   \code{X[j, ]}.
//'   Likewise, for the \code{xcorXX} functions, the \code{i}'th and
//'   \code{j}'th entry of the output is the correlation between \code{X[i, ]} 
//'   and \code{Y[j, ]}.
//' @note 
//'   \code{NA}s in \code{X} or \code{Y} will yield \code{NA}s in the correlation matrix.
//'   This also includes the diagonal unlike the behavior of 
//'   \code{stats::cor(X)}.
//' @author Anders Ellern Bilgrau <anders.ellern.bilgrau (at) gmail.com>
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix covRcpp(Rcpp::NumericMatrix & X,
                            const int norm_type = 0) {
  
  const int n = X.nrow();
  const int m = X.ncol();
  const int df = n - 1 + norm_type;

  Rcpp::NumericMatrix cov(m, m);
  
  // Degenerate case
  if (n == 0) {
    std::fill(cov.begin(), cov.end(), Rcpp::NumericVector::get_na());
    return cov; 
  }
  
  // Centering the matrix!
  X = centerNumericMatrix(X);  // Defined in aux_functions

  // Computing the covariance matrix
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
      cov(i,j) = Rcpp::sum(X(Rcpp::_, i)*X(Rcpp::_, j))/df;
      cov(j,i) = cov(i,j);
    }
  }
  
  return cov;
}


// Cross-covariance implementation in Rcpp
//' @rdname covFamily
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix xcovRcpp(Rcpp::NumericMatrix & X,
                             Rcpp::NumericMatrix & Y,
                             const int norm_type = 0) {
  
  const int n = X.nrow();
  const int m_X = X.ncol();
  const int m_Y = Y.ncol();
  const int df = n - 1 + norm_type;
  
  Rcpp::NumericMatrix cov(m_X, m_Y);
  
  // Degenerate case
  if (n == 0) {
    std::fill(cov.begin(), cov.end(), Rcpp::NumericVector::get_na());
    return cov; 
  }
  
  // Centering the matrices
  X = centerNumericMatrix(X);
  Y = centerNumericMatrix(Y);
  
  // Computing the covariance matrix
  for (int i = 0; i < m_X; ++i) {
    for (int j = 0; j < m_Y; ++j) {
      cov(i,j) = Rcpp::sum(X(Rcpp::_, i)*Y(Rcpp::_, j))/df;
    }
  }
  
  return cov;
}


// Covariance implementation in Armadillio
//' @rdname covFamily
//' @export
// [[Rcpp::export]]
arma::mat covArma(const arma::mat& X,
                  const int norm_type = 0) {
  arma::mat out(X.n_cols, X.n_cols);
  
  // Degenerate cases 
  if (X.n_cols == 0) {
    return out;
  } else if (X.n_rows == 0 || X.n_rows == 1) {
    out.fill(Rcpp::NumericVector::get_na());
  } else {
    out = arma::cov(X, norm_type);
  }
  
  return out;
}


// Cross-covariance implementation in Armadillio
//' @rdname covFamily
//' @export
// [[Rcpp::export]]
arma::mat xcovArma(const arma::mat& X,
                   const arma::mat& Y,
                   const int norm_type = 0) {
  arma::mat out(X.n_cols, Y.n_cols);
  
  // Degenerate case first
  if (X.n_cols == 0 || Y.n_cols == 0) {
    return out;
  } else if (X.n_rows == 0 || X.n_rows == 1 || Y.n_rows == 0 || Y.n_rows == 1) {
    out.fill(Rcpp::NumericVector::get_na());
  } else {
    out = arma::cov(X, Y, norm_type);
  }
  
  return out;
}

// Covariance in Eigen
//' @rdname covFamily
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd covEigen(Eigen::Map<Eigen::MatrixXd> & X,
                         const int norm_type = 0) {
  
  // Handle degenerate cases
  if (X.rows() == 0 && X.cols() > 0) {
    return Eigen::MatrixXd::Constant(X.cols(), X.cols(), 
                                     Rcpp::NumericVector::get_na());
  }
  
  // Computing degrees of freedom
  // n - 1 is the unbiased estimate whereas n is the MLE
  const int df = X.rows() - 1 + norm_type; // Subtract 1 by default
  
  // Centering the matrix
  X.rowwise() -= X.colwise().mean(); // Centering the matrix to have 0 col means
  
  // Return the X^T * X is the scatter matrix
  return X.transpose() * X / df;  // could be .adjoint()
}


// Cross-covariance in Eigen
//' @rdname covFamily
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd xcovEigen(Eigen::Map<Eigen::MatrixXd> & X,
                          Eigen::Map<Eigen::MatrixXd> & Y,
                          const int norm_type = 0) {
  
  // Handle degenerate cases
  if (X.cols() == 0 || Y.cols() == 0) {
    return Eigen::MatrixXd::Constant(0, 0, 0);
  } else if (X.rows() == 0) { // && X.cols() > 0 && Y.cols() > 0 implicit
    return Eigen::MatrixXd::Constant(X.cols(), Y.cols(),
                                     Rcpp::NumericVector::get_na());
  }
  
  // Computing degrees of freedom
  // n - 1 is the unbiased estimate whereas n is the MLE
  const int df = X.rows() - 1 + norm_type; // Subtract 1 by default
  
  // Centering matrices
  X.rowwise() -= X.colwise().mean(); // Centering the matrix to have 0 col means
  Y.rowwise() -= Y.colwise().mean();
  
  // Return the X^T * X is the scatter matrix
  return X.transpose() * Y / df;
}

