// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// corphylo_LL
double corphylo_LL(const arma::vec& par, const arma::mat& XX, const arma::mat& UU, const arma::mat& MM, const arma::mat& tau, const arma::mat& Vphy, bool REML, bool constrain_d);
RcppExport SEXP _corphyloCpp_corphylo_LL(SEXP parSEXP, SEXP XXSEXP, SEXP UUSEXP, SEXP MMSEXP, SEXP tauSEXP, SEXP VphySEXP, SEXP REMLSEXP, SEXP constrain_dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type UU(UUSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type MM(MMSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Vphy(VphySEXP);
    Rcpp::traits::input_parameter< bool >::type REML(REMLSEXP);
    Rcpp::traits::input_parameter< bool >::type constrain_d(constrain_dSEXP);
    rcpp_result_gen = Rcpp::wrap(corphylo_LL(par, XX, UU, MM, tau, Vphy, REML, constrain_d));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_corphyloCpp_corphylo_LL", (DL_FUNC) &_corphyloCpp_corphylo_LL, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_corphyloCpp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}