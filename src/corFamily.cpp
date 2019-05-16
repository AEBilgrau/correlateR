// Only include aux_functions.h which pulls in RcppArmadillo.h,
// RcppEigen.h, and Rcpp.h

#include "auxiliary_functions.h"

//// [[Rcpp::depends(RcppArmadillo)]] // Uncomment when sourceCpp()ing
//// [[Rcpp::depends(RcppEigen)]]

//' Marginal correlation matrix
//' 
//' Various workhorse functions to compute the marginal (or unconditional) 
//' correlations (and cross-correlation) estimates efficiently. 
//' They are (almost) 
//' equivalent implementations of \code{\link[stats]{cor}} in Rcpp, 
//' RcppArmadillo, and RcppEigen.
//' 
//' @rdname corFamily
//' @aliases corFamily
//'   corRcpp xcorRcpp corArma xcorArma corEigen xcorEigen
//' @param X A numeric matrix.
//' @param Y A numeric matrix of compatible dimension with the \code{X}, i.e. 
//'   \code{nrow(X)} equals \code{nrow(Y)}.
//' @return
//'   The \code{corXX} family returns a numeric correlation matrix of size 
//'   \code{ncol(X)} times \code{ncol(X)}.
//'   
//'   The \code{xcorXX} family returns a numeric cross-correlation matrix 
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
//'   \code{\link[stats]{cor}}.
//' @author Anders Ellern Bilgrau <anders.ellern.bilgrau (at) gmail.com>
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix corRcpp(Rcpp::NumericMatrix & X) {
  
  const int m = X.ncol();
  const int n = X.nrow();
  
  // Centering the matrix
  X = centerNumericMatrix(X);
  
  Rcpp::NumericMatrix cor(m, m);
  
  // Degenerate case
  if (n == 0) {
    std::fill(cor.begin(), cor.end(), Rcpp::NumericVector::get_na());
    return cor; 
  }
  
  // Compute 1 over the sample standard deviation
  Rcpp::NumericVector inv_sqrt_ss(m);
  for (int i = 0; i < m; ++i) {
    inv_sqrt_ss(i) = 1/sqrt(Rcpp::sum(X(Rcpp::_, i)*X(Rcpp::_, i)));
  }
  
  // Computing the correlation matrix
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
      cor(i, j) = Rcpp::sum(X(Rcpp::_,i)*X(Rcpp::_,j)) *
        inv_sqrt_ss(i) * inv_sqrt_ss(j);
      cor(j, i) = cor(i, j);
    }
  }

  return cor;
}


// Cross-correlation implementation in Rcpp
//' @rdname corFamily
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix xcorRcpp(Rcpp::NumericMatrix & X,
                             Rcpp::NumericMatrix & Y) {
  
  const int m_X = X.ncol();
  const int m_Y = Y.ncol();
  const int n = X.nrow();
  
  // Centering the matrices
  X = centerNumericMatrix(X);
  Y = centerNumericMatrix(Y);
  
  Rcpp::NumericMatrix cor(m_X, m_Y);
  
  // Degenerate case
  if (n == 0) {
    std::fill(cor.begin(), cor.end(), Rcpp::NumericVector::get_na());
    return cor; 
  }
  
  // Compute 1 over square root the sum of squares
  Rcpp::NumericVector inv_sqrt_ss_X(m_X);
  for (int i = 0; i < m_X; ++i) {
    inv_sqrt_ss_X(i) = 1/sqrt(Rcpp::sum(X(Rcpp::_, i)*X(Rcpp::_, i)));
  }
  Rcpp::NumericVector inv_sqrt_ss_Y(m_Y);
  for (int i = 0; i < m_Y; ++i) {
    inv_sqrt_ss_Y(i) = 1/sqrt(Rcpp::sum(Y(Rcpp::_, i)*Y(Rcpp::_, i)));
  }
  
  // Computing the cross-correlation matrix
  for (int i = 0; i < m_X; ++i) {
    for (int j = 0; j < m_Y; ++j) {
      cor(i, j) = Rcpp::sum(X(Rcpp::_, i)*Y(Rcpp::_, j)) *
        inv_sqrt_ss_X(i) * inv_sqrt_ss_Y(j);
    }
  }

  return cor;
}



// Correlation implementation in armadillo
//' @rdname corFamily
//' @export
// [[Rcpp::export]]
arma::mat corArma(const arma::mat & X) {
  arma::mat out(X.n_cols, X.n_cols);
  
  // Degenerate cases 
  if (X.n_cols == 0) {
    return out;
  } else if (X.n_rows == 0 || X.n_rows == 1) {
    out.fill(Rcpp::NumericVector::get_na());
  } else {
    out = arma::cor(X, 0);
  }

  out.diag() /= out.diag();   // Ensure 1 in the diagonal (if not NaN)
  return out;
}

// Cross-correlation implementation in armadillo
//' @rdname corFamily
//' @export
// [[Rcpp::export]]
arma::mat xcorArma(const arma::mat& X,
                   const arma::mat& Y) {
  arma::mat out(X.n_cols, Y.n_cols);
  
  // Degenerate case first
  if (X.n_cols == 0 || Y.n_cols == 0) {
    return out;
  } else if (X.n_rows == 0 || X.n_rows == 1 || Y.n_rows == 0 || Y.n_rows == 1) {
    out.fill(Rcpp::NumericVector::get_na());
  } else {
    out = arma::cor(X, Y, 0);
  }
  
  return out;
}




// Correlation implementation in Eigen
//' @rdname corFamily
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd corEigen(Eigen::Map<Eigen::MatrixXd> & X) {
  
  // Handle degenerate cases
  if (X.rows() == 0 && X.cols() > 0) {
    return Eigen::MatrixXd::Constant(X.cols(), X.cols(), 
                                     Rcpp::NumericVector::get_na());
  }
  
  // Computing degrees of freedom
  // n - 1 is the unbiased estimate whereas n is the MLE
  const int df = X.rows() - 1; // Subtract 1 by default
  
  X.rowwise() -= X.colwise().mean();  // Centering
  
  Eigen::MatrixXd cor = X.transpose() * X / df;   // The covariance matrix
  
  // Get 1 over the standard deviations
  Eigen::VectorXd inv_sds = cor.diagonal().array().sqrt().inverse();
  
  // Scale the covariance matrix
  cor = cor.cwiseProduct(inv_sds * inv_sds.transpose());
  
  return cor;
}

// Cross-correlation implementation in Eigen
//' @rdname corFamily
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd xcorEigen(Eigen::Map<Eigen::MatrixXd> & X,
                          Eigen::Map<Eigen::MatrixXd> & Y) {
  
  // Handle degenerate cases
  if (X.cols() == 0 || Y.cols() == 0) {
    return Eigen::MatrixXd::Constant(0, 0, 0);
  } else if (X.rows() == 0) { // && X.cols() > 0 && Y.cols() > 0 implicit
    return Eigen::MatrixXd::Constant(X.cols(), Y.cols(),
                                     Rcpp::NumericVector::get_na());
  }
  
  // Computing degrees of freedom
  // n - 1 is the unbiased estimate whereas n is the MLE
  const int df = X.rows() - 1; // Subtract 1 by default
  
  // Centering matrices
  Y.rowwise() -= Y.colwise().mean();
  X.rowwise() -= X.colwise().mean();
  
  // The covariance matrix
  Eigen::MatrixXd cor = X.transpose() * Y / df;
  
  // Compute 1 over the standard deviations of X and Y
  Eigen::VectorXd inv_sds_X = (X.colwise().norm()/sqrt(df)).array().inverse();
  Eigen::VectorXd inv_sds_Y = (Y.colwise().norm()/sqrt(df)).array().inverse();
  
  // Scale the covariance matrix
  cor = cor.cwiseProduct(inv_sds_X * inv_sds_Y.transpose());
  return cor;
}


