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
//' equivalent implementations of \code{\link{stats::cov}} in Rcpp, 
//' RcppArmadillo, and RcppEigen.
//' 
//' @rdname covFamily
//' @aliases covFamily
//'   corRcpp xcorRcpp corArma xcorArma corEigen xcorEigen
//' @param X A numeric matrix.
//' @param Y A numeric matrix of compatible dimension with the \code{X}, i.e. 
//'   \code{nrow(X)} equals \code{nrow(Y)}.
//' @param norm.type an integer of length one giving the estimator. The 
//'   default \code{0L} gives the unbised estimate while \code{1L} gives the 
//'   MLE.
//' @return
//'   The \code{corXX} familiy returns a numeric correlation matrix of size 
//'   \code{ncol(X)} times \code{ncol(X)}.
//'   
//'   The \code{xcorXX} family returns a numeric cross-covariance matrix 
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
//'   NA's in \code{X} or \code{Y} will yield NA's in the correlation matrix.
//'   This also includes the diagonal unlike the behaviour of 
//'   \code{stats::cor(X)}.
//' @author Anders Ellern Bilgrau <abilgrau (at) math.aau.dk>
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix covRcpp(Rcpp::NumericMatrix & X,
                            const int norm_type = 0) {
  
  const int n = X.nrow();
  const int m = X.ncol();
  const int df = n - 1 + norm_type;
  
  // Centering the matrix!
  X = centerNumericMatrix(X);  // Defined in aux_functions

  // Computing the covariance matrix
  Rcpp::NumericMatrix cov(m, m);
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
  
  // Centering the matrices
  X = centerNumericMatrix(X);
  Y = centerNumericMatrix(Y);
  
  // Computing the covariance matrix
  Rcpp::NumericMatrix cov(m_X, m_Y);
  for (int i = 0; i < m_X; ++i) {
    for (int j = 0; j < m_Y; ++j) {
      cov(i,j) = Rcpp::sum(X(Rcpp::_, i)*Y(Rcpp::_, j))/df;
    }
  }
  
  return cov;
}


// Covariance "implementation"" in Armadillio
//' @rdname covFamily
//' @export
// [[Rcpp::export]]
arma::mat covArma(const arma::mat & X,
                  const int norm_type = 0) {
  return arma::cov(X, norm_type);
}


// Cross-covariance "implementation"" in Armadillio
//' @rdname covFamily
//' @export
// [[Rcpp::export]]
arma::mat xcovArma(const arma::mat & X,
                       const arma::mat & Y,
                       const int norm_type = 0) {
  return arma::cov(X, Y, norm_type);
}


// covariance in Eigen
//' @rdname covFamily
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd covEigen(Eigen::Map<Eigen::MatrixXd> & X,
                         const int norm_type = 0) {

    // Computing degrees of freedom
    // n - 1 is the unbiased estimate whereas n is the MLE
    const int df = X.rows() - 1 + norm_type; // Subtract 1 by default

    // Centering the matrix
    X.rowwise() -= X.colwise().mean();  

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

    // Computing degrees of freedom
    // n - 1 is the unbiased estimate whereas n is the MLE
    const int df = X.rows() - 1 + norm_type; // Subtract 1 by default

    // Centering matrices
    Y.rowwise() -= Y.colwise().mean();
    X.rowwise() -= X.colwise().mean();

    // Return the X^T * X is the scatter matrix
    return X.transpose() * Y / df;
}



/*** R
library("microbenchmark")

X <- replicate(10, rnorm(50))
dimnames(X) <- list(paste0("obs", 1:nrow(X)), paste0("dim", 1:ncol(X)))
Y <- replicate(15, rnorm(50))
dimnames(Y) <- list(paste0("obs", 1:nrow(Y)), paste0("dim", 1:ncol(Y)))

covArma(X, 0)
covRcpp(X, 0)
covEigen(X, 0)
all.equal(covArma(X, 0), covRcpp(X, 0))
all.equal(covArma(X, 1), covRcpp(X, 1))

all(crosscovRcpp(X, Y, 1) == crosscovArma(X, Y, 1))
print(microbenchmark(crosscovArma(X, Y, 0), 
                     crosscovEigen(X, Y, 0),
                     crosscovRcpp(X, Y, 0), times = 10), order = "median")
*/

