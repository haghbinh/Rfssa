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
// HLinprod
double HLinprod(arma::mat x, arma::mat y, arma::mat G);
RcppExport SEXP _Rfssa_HLinprod(SEXP xSEXP, SEXP ySEXP, SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(HLinprod(x, y, G));
    return rcpp_result_gen;
END_RCPP
}
// HpLinprod
double HpLinprod(std::vector<arma::mat> X, std::vector<arma::mat> Y, std::vector<arma::mat> G, int p);
RcppExport SEXP _Rfssa_HpLinprod(SEXP XSEXP, SEXP YSEXP, SEXP GSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type X(XSEXP);
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(HpLinprod(X, Y, G, p));
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
// SSM
NumericMatrix SSM(int K, int L, int d_tilde, int p, std::vector<NumericMatrix> B, NumericMatrix shifter);
RcppExport SEXP _Rfssa_SSM(SEXP KSEXP, SEXP LSEXP, SEXP d_tildeSEXP, SEXP pSEXP, SEXP BSEXP, SEXP shifterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type d_tilde(d_tildeSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< std::vector<NumericMatrix> >::type B(BSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type shifter(shifterSEXP);
    rcpp_result_gen = Rcpp::wrap(SSM(K, L, d_tilde, p, B, shifter));
    return rcpp_result_gen;
END_RCPP
}
// Gramm
NumericMatrix Gramm(int K, int L, int p, int d_tilde, std::vector<NumericMatrix> A, NumericMatrix shifter, NumericMatrix d);
RcppExport SEXP _Rfssa_Gramm(SEXP KSEXP, SEXP LSEXP, SEXP pSEXP, SEXP d_tildeSEXP, SEXP ASEXP, SEXP shifterSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type d_tilde(d_tildeSEXP);
    Rcpp::traits::input_parameter< std::vector<NumericMatrix> >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type shifter(shifterSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(Gramm(K, L, p, d_tilde, A, shifter, d));
    return rcpp_result_gen;
END_RCPP
}
// mwinprod
double mwinprod(std::vector<arma::mat> X, std::vector<arma::mat> Y, arma::vec w, std::vector<arma::mat> G, int p);
RcppExport SEXP _Rfssa_mwinprod(SEXP XSEXP, SEXP YSEXP, SEXP wSEXP, SEXP GSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type X(XSEXP);
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(mwinprod(X, Y, w, G, p));
    return rcpp_result_gen;
END_RCPP
}
// winprod
double winprod(arma::mat x, arma::mat y, arma::vec w, arma::mat G);
RcppExport SEXP _Rfssa_winprod(SEXP xSEXP, SEXP ySEXP, SEXP wSEXP, SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(winprod(x, y, w, G));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Rfssa_H", (DL_FUNC) &_Rfssa_H, 1},
    {"_Rfssa_HLinprod", (DL_FUNC) &_Rfssa_HLinprod, 3},
    {"_Rfssa_HpLinprod", (DL_FUNC) &_Rfssa_HpLinprod, 4},
    {"_Rfssa_Cofmat", (DL_FUNC) &_Rfssa_Cofmat, 3},
    {"_Rfssa_SS", (DL_FUNC) &_Rfssa_SS, 4},
    {"_Rfssa_SSM", (DL_FUNC) &_Rfssa_SSM, 6},
    {"_Rfssa_Gramm", (DL_FUNC) &_Rfssa_Gramm, 7},
    {"_Rfssa_mwinprod", (DL_FUNC) &_Rfssa_mwinprod, 5},
    {"_Rfssa_winprod", (DL_FUNC) &_Rfssa_winprod, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_Rfssa(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
