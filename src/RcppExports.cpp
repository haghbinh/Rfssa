// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// H
arma::mat H(arma::mat A);
RcppExport SEXP _Rfssa_H(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(H(A));
    return rcpp_result_gen;
END_RCPP
}
// mod
NumericVector mod(int K, int W);
RcppExport SEXP _Rfssa_mod(SEXP KSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(mod(K, W));
    return rcpp_result_gen;
END_RCPP
}
// Csij
double Csij(int i, int j, int K, int L, NumericMatrix B);
RcppExport SEXP _Rfssa_Csij(SEXP iSEXP, SEXP jSEXP, SEXP KSEXP, SEXP LSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(Csij(i, j, K, L, B));
    return rcpp_result_gen;
END_RCPP
}
// SS
NumericMatrix SS(int K, int L, NumericMatrix B, int d);
RcppExport SEXP _Rfssa_SS(SEXP KSEXP, SEXP LSEXP, SEXP BSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(SS(K, L, B, d));
    return rcpp_result_gen;
END_RCPP
}
// Cofmat
NumericMatrix Cofmat(int d, int L, NumericVector cx);
RcppExport SEXP _Rfssa_Cofmat(SEXP dSEXP, SEXP LSEXP, SEXP cxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cx(cxSEXP);
    rcpp_result_gen = Rcpp::wrap(Cofmat(d, L, cx));
    return rcpp_result_gen;
END_RCPP
}
// Gram
NumericMatrix Gram(int K, int L, NumericMatrix A, int d);
RcppExport SEXP _Rfssa_Gram(SEXP KSEXP, SEXP LSEXP, SEXP ASEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(Gram(K, L, A, d));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Rfssa_H", (DL_FUNC) &_Rfssa_H, 1},
    {"_Rfssa_mod", (DL_FUNC) &_Rfssa_mod, 2},
    {"_Rfssa_Csij", (DL_FUNC) &_Rfssa_Csij, 5},
    {"_Rfssa_SS", (DL_FUNC) &_Rfssa_SS, 4},
    {"_Rfssa_Cofmat", (DL_FUNC) &_Rfssa_Cofmat, 3},
    {"_Rfssa_Gram", (DL_FUNC) &_Rfssa_Gram, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_Rfssa(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
