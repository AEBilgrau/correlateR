// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// rcm_logdetPsiPlusS_arma
Rcpp::NumericVector rcm_logdetPsiPlusS_arma(const arma::mat& Psi, const Rcpp::List& S_list);
RcppExport SEXP _correlateR_rcm_logdetPsiPlusS_arma(SEXP PsiSEXP, SEXP S_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type S_list(S_listSEXP);
    rcpp_result_gen = Rcpp::wrap(rcm_logdetPsiPlusS_arma(Psi, S_list));
    return rcpp_result_gen;
END_RCPP
}
// rcm_loglik_arma
double rcm_loglik_arma(const arma::mat& Psi, const double nu, const Rcpp::List& S_list, const Rcpp::NumericVector& ns);
RcppExport SEXP _correlateR_rcm_loglik_arma(SEXP PsiSEXP, SEXP nuSEXP, SEXP S_listSEXP, SEXP nsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< const double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type S_list(S_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type ns(nsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcm_loglik_arma(Psi, nu, S_list, ns));
    return rcpp_result_gen;
END_RCPP
}
// rcm_loglik_nu_arma
double rcm_loglik_nu_arma(const double& logdetPsi, const double nu, const Rcpp::NumericVector& logdetPsiPlusS, const Rcpp::NumericVector& ns, const int p);
RcppExport SEXP _correlateR_rcm_loglik_nu_arma(SEXP logdetPsiSEXP, SEXP nuSEXP, SEXP logdetPsiPlusSSEXP, SEXP nsSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type logdetPsi(logdetPsiSEXP);
    Rcpp::traits::input_parameter< const double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type logdetPsiPlusS(logdetPsiPlusSSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(rcm_loglik_nu_arma(logdetPsi, nu, logdetPsiPlusS, ns, p));
    return rcpp_result_gen;
END_RCPP
}
// rcm_em_step_arma
arma::mat rcm_em_step_arma(const arma::mat& Psi, const double nu, const Rcpp::List& S_list, const Rcpp::NumericVector& ns);
RcppExport SEXP _correlateR_rcm_em_step_arma(SEXP PsiSEXP, SEXP nuSEXP, SEXP S_listSEXP, SEXP nsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< const double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type S_list(S_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type ns(nsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcm_em_step_arma(Psi, nu, S_list, ns));
    return rcpp_result_gen;
END_RCPP
}
// lgammap
Rcpp::NumericVector lgammap(const Rcpp::NumericVector& x, const int p);
RcppExport SEXP _correlateR_lgammap(SEXP xSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(lgammap(x, p));
    return rcpp_result_gen;
END_RCPP
}
// logdet_arma
arma::vec logdet_arma(const arma::mat& x);
RcppExport SEXP _correlateR_logdet_arma(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logdet_arma(x));
    return rcpp_result_gen;
END_RCPP
}
// corRcpp
Rcpp::NumericMatrix corRcpp(Rcpp::NumericMatrix& X);
RcppExport SEXP _correlateR_corRcpp(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(corRcpp(X));
    return rcpp_result_gen;
END_RCPP
}
// xcorRcpp
Rcpp::NumericMatrix xcorRcpp(Rcpp::NumericMatrix& X, Rcpp::NumericMatrix& Y);
RcppExport SEXP _correlateR_xcorRcpp(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(xcorRcpp(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// corArma
arma::mat corArma(const arma::mat& X);
RcppExport SEXP _correlateR_corArma(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(corArma(X));
    return rcpp_result_gen;
END_RCPP
}
// xcorArma
arma::mat xcorArma(const arma::mat& X, const arma::mat& Y);
RcppExport SEXP _correlateR_xcorArma(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(xcorArma(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// corEigen
Eigen::MatrixXd corEigen(Eigen::Map<Eigen::MatrixXd>& X);
RcppExport SEXP _correlateR_corEigen(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(corEigen(X));
    return rcpp_result_gen;
END_RCPP
}
// xcorEigen
Eigen::MatrixXd xcorEigen(Eigen::Map<Eigen::MatrixXd>& X, Eigen::Map<Eigen::MatrixXd>& Y);
RcppExport SEXP _correlateR_xcorEigen(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(xcorEigen(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// cov2corArma
arma::mat cov2corArma(arma::mat S);
RcppExport SEXP _correlateR_cov2corArma(SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(cov2corArma(S));
    return rcpp_result_gen;
END_RCPP
}
// covRcpp
Rcpp::NumericMatrix covRcpp(Rcpp::NumericMatrix& X, const int norm_type);
RcppExport SEXP _correlateR_covRcpp(SEXP XSEXP, SEXP norm_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int >::type norm_type(norm_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(covRcpp(X, norm_type));
    return rcpp_result_gen;
END_RCPP
}
// xcovRcpp
Rcpp::NumericMatrix xcovRcpp(Rcpp::NumericMatrix& X, Rcpp::NumericMatrix& Y, const int norm_type);
RcppExport SEXP _correlateR_xcovRcpp(SEXP XSEXP, SEXP YSEXP, SEXP norm_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const int >::type norm_type(norm_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(xcovRcpp(X, Y, norm_type));
    return rcpp_result_gen;
END_RCPP
}
// covArma
arma::mat covArma(const arma::mat& X, const int norm_type);
RcppExport SEXP _correlateR_covArma(SEXP XSEXP, SEXP norm_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int >::type norm_type(norm_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(covArma(X, norm_type));
    return rcpp_result_gen;
END_RCPP
}
// xcovArma
arma::mat xcovArma(const arma::mat& X, const arma::mat& Y, const int norm_type);
RcppExport SEXP _correlateR_xcovArma(SEXP XSEXP, SEXP YSEXP, SEXP norm_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const int >::type norm_type(norm_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(xcovArma(X, Y, norm_type));
    return rcpp_result_gen;
END_RCPP
}
// covEigen
Eigen::MatrixXd covEigen(Eigen::Map<Eigen::MatrixXd>& X, const int norm_type);
RcppExport SEXP _correlateR_covEigen(SEXP XSEXP, SEXP norm_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int >::type norm_type(norm_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(covEigen(X, norm_type));
    return rcpp_result_gen;
END_RCPP
}
// xcovEigen
Eigen::MatrixXd xcovEigen(Eigen::Map<Eigen::MatrixXd>& X, Eigen::Map<Eigen::MatrixXd>& Y, const int norm_type);
RcppExport SEXP _correlateR_xcovEigen(SEXP XSEXP, SEXP YSEXP, SEXP norm_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const int >::type norm_type(norm_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(xcovEigen(X, Y, norm_type));
    return rcpp_result_gen;
END_RCPP
}
// pcorArma
arma::mat pcorArma(const arma::mat& X, const arma::uvec& z);
RcppExport SEXP _correlateR_pcorArma(SEXP XSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(pcorArma(X, z));
    return rcpp_result_gen;
END_RCPP
}
// pcovArma
arma::mat pcovArma(const arma::mat& X, const arma::uvec& z);
RcppExport SEXP _correlateR_pcovArma(SEXP XSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(pcovArma(X, z));
    return rcpp_result_gen;
END_RCPP
}
// pool
arma::mat pool(const Rcpp::List& S_list, const Rcpp::NumericVector ns, const int norm_type);
RcppExport SEXP _correlateR_pool(SEXP S_listSEXP, SEXP nsSEXP, SEXP norm_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type S_list(S_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const int >::type norm_type(norm_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(pool(S_list, ns, norm_type));
    return rcpp_result_gen;
END_RCPP
}
// pxcorArma
arma::mat pxcorArma(const arma::mat& X, const arma::mat& Y, const arma::mat& Z);
RcppExport SEXP _correlateR_pxcorArma(SEXP XSEXP, SEXP YSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(pxcorArma(X, Y, Z));
    return rcpp_result_gen;
END_RCPP
}
// pxcovArma
arma::mat pxcovArma(const arma::mat& X, const arma::mat& Y, const arma::mat& Z, const int norm_type);
RcppExport SEXP _correlateR_pxcovArma(SEXP XSEXP, SEXP YSEXP, SEXP ZSEXP, SEXP norm_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const int >::type norm_type(norm_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(pxcovArma(X, Y, Z, norm_type));
    return rcpp_result_gen;
END_RCPP
}
// rwishartArma
arma::cube rwishartArma(const int n, const arma::mat& sigma, const double nu);
RcppExport SEXP _correlateR_rwishartArma(SEXP nSEXP, SEXP sigmaSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(rwishartArma(n, sigma, nu));
    return rcpp_result_gen;
END_RCPP
}
// rinvwishartArma
arma::cube rinvwishartArma(const int n, const arma::mat& psi, const double nu);
RcppExport SEXP _correlateR_rinvwishartArma(SEXP nSEXP, SEXP psiSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(rinvwishartArma(n, psi, nu));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_correlateR_rcm_logdetPsiPlusS_arma", (DL_FUNC) &_correlateR_rcm_logdetPsiPlusS_arma, 2},
    {"_correlateR_rcm_loglik_arma", (DL_FUNC) &_correlateR_rcm_loglik_arma, 4},
    {"_correlateR_rcm_loglik_nu_arma", (DL_FUNC) &_correlateR_rcm_loglik_nu_arma, 5},
    {"_correlateR_rcm_em_step_arma", (DL_FUNC) &_correlateR_rcm_em_step_arma, 4},
    {"_correlateR_lgammap", (DL_FUNC) &_correlateR_lgammap, 2},
    {"_correlateR_logdet_arma", (DL_FUNC) &_correlateR_logdet_arma, 1},
    {"_correlateR_corRcpp", (DL_FUNC) &_correlateR_corRcpp, 1},
    {"_correlateR_xcorRcpp", (DL_FUNC) &_correlateR_xcorRcpp, 2},
    {"_correlateR_corArma", (DL_FUNC) &_correlateR_corArma, 1},
    {"_correlateR_xcorArma", (DL_FUNC) &_correlateR_xcorArma, 2},
    {"_correlateR_corEigen", (DL_FUNC) &_correlateR_corEigen, 1},
    {"_correlateR_xcorEigen", (DL_FUNC) &_correlateR_xcorEigen, 2},
    {"_correlateR_cov2corArma", (DL_FUNC) &_correlateR_cov2corArma, 1},
    {"_correlateR_covRcpp", (DL_FUNC) &_correlateR_covRcpp, 2},
    {"_correlateR_xcovRcpp", (DL_FUNC) &_correlateR_xcovRcpp, 3},
    {"_correlateR_covArma", (DL_FUNC) &_correlateR_covArma, 2},
    {"_correlateR_xcovArma", (DL_FUNC) &_correlateR_xcovArma, 3},
    {"_correlateR_covEigen", (DL_FUNC) &_correlateR_covEigen, 2},
    {"_correlateR_xcovEigen", (DL_FUNC) &_correlateR_xcovEigen, 3},
    {"_correlateR_pcorArma", (DL_FUNC) &_correlateR_pcorArma, 2},
    {"_correlateR_pcovArma", (DL_FUNC) &_correlateR_pcovArma, 2},
    {"_correlateR_pool", (DL_FUNC) &_correlateR_pool, 3},
    {"_correlateR_pxcorArma", (DL_FUNC) &_correlateR_pxcorArma, 3},
    {"_correlateR_pxcovArma", (DL_FUNC) &_correlateR_pxcovArma, 4},
    {"_correlateR_rwishartArma", (DL_FUNC) &_correlateR_rwishartArma, 3},
    {"_correlateR_rinvwishartArma", (DL_FUNC) &_correlateR_rinvwishartArma, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_correlateR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
