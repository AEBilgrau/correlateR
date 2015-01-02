// Only include aux_functions.h which pulls in RcppArmadillo.h,
// RcppEigen.h, and Rcpp.h

#include "auxiliary_functions.h"

//' The RCM log-likelihood function
//' 
//' Fast evaluation of the RCM log-likelihood function.
//' 
//' @param Psi A \code{p} times \code{p} numeric positive semi-definte matrix.
//' @param nu A numeric of length 1 giving the degrees of freedom.
//' @param S_list A \code{list} of scatter matrices of the same size as 
//'   \code{Psi} for each group.
//' @param ns A numeric of the same length as \code{S_list} giving the 
//'   number of samples in each group.
//' @return The value of the log-likelihood.
//' @author Anders Ellern Bilgrau <abilgrau (at) math.aau.dk>
//' @examples
//' ns <-  c(5, 5, 5)
//' S <- createRCMData(ns = ns, psi = diag(4), nu = 30)
//' correlateR:::rcm_loglik_arma(Psi = diag(4), nu = 15, S_list = S, ns = ns)
//' @keywords internal
// [[Rcpp::export]]
double rcm_loglik_arma(const arma::mat & Psi, 
                       const double nu, 
                       const Rcpp::List & S_list,  
                       const Rcpp::NumericVector & ns) { 
  const int k = S_list.size();
  const int p = Psi.n_rows;
  const Rcpp::NumericVector nu_half(1, nu/2.0f); // Vector of length 1
  const Rcpp::NumericVector cs = (nu + ns)/2.0f;
  const double logdetPsi = logdet_arma( Psi )(0);
  
  arma::vec logdetPsiPlusS(k);
  for (int i = 0; i < k; ++i) {
    arma::mat Si = S_list[i];
    logdetPsiPlusS(i) = logdet_arma( Psi + Si )(0);
  }
  const double c1 = sum(ns * p)/2.0f * log(arma::datum::pi);
  const double t1 = k*nu_half(0)*logdetPsi;
  const double t2 = sum( lgammap(cs, p) );
  const double t3 = -sum(cs*Rcpp::as<Rcpp::NumericVector>(
    Rcpp::wrap(logdetPsiPlusS)));
  const double t4 = -k * lgammap(nu_half, p)(0);
  return c1 + t1 + t2 + t3 + t4;
}







//' The RCM EM-step
//' 
//' A armadillo-based function to perform the E and M step in the 
//' EM algorithm of the RCM. This functions assumes \code{nu} to be fixed.
//' 
//' @param Psi A numeric matrix.
//' @param nu A numeric of length 1 giving the degrees of freedom in the RCM.
//' @param S_list A \code{list} of scatter matrices for each dataset/group
//'   of the same size a \code{Psi}.
//' @param ns A numeric vector the same lengths as \code{S_list} giving the
//'   number of samples for each dataset.
//' @return A numeric matrix the same size as \code{Psi} giving the updated
//'   \code{Psi}.
//' @examples
//' ns <-  c(5, 5, 5)
//' S <- createRCMData(ns = ns, psi = diag(4), nu = 30)
//' correlateR:::rcm_em_step_arma(Psi = diag(4), nu = 15, S_list = S, ns = ns)
//' @keywords internal 
// [[Rcpp::export]]
arma::mat rcm_em_step_arma(const arma::mat & Psi, 
                           const double nu, 
                           const Rcpp::List & S_list,  
                           const Rcpp::NumericVector & ns) {
  const int k = S_list.size();
  const int p = Psi.n_rows;
  const double c = k*nu;
  arma::mat inv_ans(p, p, arma::fill::zeros);
  for (int i = 0; i < k; ++i) {
      double fac = ns[i] + nu;
      arma::mat Si = S_list[i];
      inv_ans += fac * arma::inv( Psi + Si );
  }
  return c * arma::inv(inv_ans);
}






/*** R

set.seed(1)
ns <-  c(5, 5, 5)
S <- createRCMData(ns = ns, psi = diag(4), nu = 30)
correlateR:::rcm_loglik_arma(Psi = diag(4), nu = 15, S_list = S, ns = ns)
 
rm(list = ls())

library("correlateR")
ns <-  c(200, 200)  # Study sizes
Psi <- diag(3)  # True Psi
nu <- 1000  # inter study homogenetiy
p <- nrow(Psi)
S <- createRCMData(ns = ns, psi = Psi, nu = nu) # Observed scatter matrices

L <- function(nu) correlateR:::rcm_loglik_arma(Psi = Psi, nu = nu, S_list = S, ns = ns)
nu.est <- nu.est2 <- c()
for (j in 0:200) {
  S <- createRCMData(ns = ns, psi = Psi, nu = nu) # Observed scatter matrices
  nus <- seq(p + 1 + .001, 1e3, l = 200)
  eva <- eva2 <- c()
  for (i in seq_along(nus)) {
    eva[i] <- L(nus[i])
  }  
  par(new = as.logical(j))
  plot(nus, eva, type = "l", axes = FALSE, xlim = c(0, max(nus)), 
       col = "#00000050", xlab = "", ylab = "")
  nu.est[j+1] <- correlateR:::rcm_get_nu(Psi = Psi, S_list = S, ns = ns)$maximum
  nu.est2[j+1] <- correlateR:::rcm_get_nu2(Psi = Psi, S_list = S, ns = ns)$par
}
axis(1)
rug(nu.est, col = "red", lwd = 2, ticksize = 0.1)
rug(nu.est, col = "red", lwd = 2, side = 3, ticksize = 0.01)
rug(nu.est2, col = "blue", lwd = 2, ticksize = 0.06)
abline(v = nu)
abline(v = mean(nu.est), col = "red")
abline(v = mean(nu.est2), col = "blue")


hist(nu.est, breaks = 50)


rm(list = ls())

library("correlateR")
ns <-  c(200, 200)  # Study sizes
Psi <- diag(3)  # True Psi
nu <- 1000  # inter study homogenetiy
p <- nrow(Psi)
S <- createRCMData(ns = ns, psi = Psi, nu = nu) # Observed scatter matrices


PsiHat <- diag(3)#correlateR:::pool(S, ns)
a <- correlateR:::rcm_loglik_arma(Psi = PsiHat, ns = ns, S_list = S, nu = nu) 
cat(sprintf("L = %0.10f\tdelta = %0.10f\n", b, NA))
for (i in 1:20) {
  PsiHat <- correlateR:::rcm_em_step_arma(Psi = PsiHat, nu = nu, S_list = S, ns = ns)
  b <- correlateR:::rcm_loglik_arma(Psi = PsiHat, ns = ns, S_list = S, nu = nu)
  cat(sprintf("L = %0.10f\tdelta = %0.10f\n", b, b-a))
  a <- b
}





emStep <- function(Psi, nu, S, ns) {
  p <- nrow(Psi)
  k <- length(S)
  co <- 1/(k*nu)
  tmp <- matrix(0, p, p)
  for (i in 1:k) {
    tmp <- tmp + (ns[i]+nu)*solve(1/(nu - p - 1)*Psi + S[[i]] )    
  }
  return(solve(co*tmp))
}

correlateR:::rcm_loglik_arma(Psi = Psi, ns = ns, S_list = S, nu = nu) 
correlateR:::rcm_loglik_arma(Psi = correlateR:::pool(S, ns), ns = ns, S_list = S, nu = nu) 

nu.new <- correlateR:::rcm_get_nu2(Psi = Psi, S_list = S, ns = ns)$estimate
correlateR:::rcm_loglik_arma(Psi = Psi, ns = ns, S_list = S, nu = nu) 
correlateR:::rcm_loglik_arma(Psi = Psi, ns = ns, S_list = S, nu = nu.new) 


(Psi3 <- emStep(Psi = Psi, nu = nu, S = S, ns = ns))
correlateR:::rcm_loglik_arma(Psi = Psi, ns = ns, S_list = S, nu = nu) 
correlateR:::rcm_loglik_arma(Psi = Psi3, ns = ns, S_list = S, nu = nu)



emStep(Psi, nu, S, ns)
correlateR:::rcm_em_step_arma(Psi = Psi, nu = nu, S_list = S, ns = ns)


*/

