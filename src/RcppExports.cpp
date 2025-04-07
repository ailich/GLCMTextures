// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// C_make_glcm_counts
IntegerMatrix C_make_glcm_counts(IntegerMatrix x, int n_levels, IntegerVector shift, bool na_rm);
RcppExport SEXP _GLCMTextures_C_make_glcm_counts(SEXP xSEXP, SEXP n_levelsSEXP, SEXP shiftSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n_levels(n_levelsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type shift(shiftSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(C_make_glcm_counts(x, n_levels, shift, na_rm));
    return rcpp_result_gen;
END_RCPP
}
// C_make_glcm
arma::mat C_make_glcm(IntegerMatrix x, int n_levels, IntegerVector shift, bool na_rm);
RcppExport SEXP _GLCMTextures_C_make_glcm(SEXP xSEXP, SEXP n_levelsSEXP, SEXP shiftSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n_levels(n_levelsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type shift(shiftSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(C_make_glcm(x, n_levels, shift, na_rm));
    return rcpp_result_gen;
END_RCPP
}
// C_GLSV
NumericVector C_GLSV(arma::mat Pij, int n_levels);
RcppExport SEXP _GLCMTextures_C_GLSV(SEXP PijSEXP, SEXP n_levelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Pij(PijSEXP);
    Rcpp::traits::input_parameter< int >::type n_levels(n_levelsSEXP);
    rcpp_result_gen = Rcpp::wrap(C_GLSV(Pij, n_levels));
    return rcpp_result_gen;
END_RCPP
}
// C_glcm_metrics
NumericVector C_glcm_metrics(arma::mat Pij, arma::mat i_mat, arma::mat j_mat, int n_levels, NumericVector k_vals, CharacterVector metrics, bool impute_corr);
RcppExport SEXP _GLCMTextures_C_glcm_metrics(SEXP PijSEXP, SEXP i_matSEXP, SEXP j_matSEXP, SEXP n_levelsSEXP, SEXP k_valsSEXP, SEXP metricsSEXP, SEXP impute_corrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Pij(PijSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type i_mat(i_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type j_mat(j_matSEXP);
    Rcpp::traits::input_parameter< int >::type n_levels(n_levelsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type k_vals(k_valsSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type metrics(metricsSEXP);
    Rcpp::traits::input_parameter< bool >::type impute_corr(impute_corrSEXP);
    rcpp_result_gen = Rcpp::wrap(C_glcm_metrics(Pij, i_mat, j_mat, n_levels, k_vals, metrics, impute_corr));
    return rcpp_result_gen;
END_RCPP
}
// C_glcm_textures_helper
NumericMatrix C_glcm_textures_helper(IntegerVector x, IntegerVector w2, int n_levels, IntegerVector shift, CharacterVector metrics, bool na_rm, bool impute_corr, size_t ni, size_t nw);
RcppExport SEXP _GLCMTextures_C_glcm_textures_helper(SEXP xSEXP, SEXP w2SEXP, SEXP n_levelsSEXP, SEXP shiftSEXP, SEXP metricsSEXP, SEXP na_rmSEXP, SEXP impute_corrSEXP, SEXP niSEXP, SEXP nwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type w2(w2SEXP);
    Rcpp::traits::input_parameter< int >::type n_levels(n_levelsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type shift(shiftSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type metrics(metricsSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< bool >::type impute_corr(impute_corrSEXP);
    Rcpp::traits::input_parameter< size_t >::type ni(niSEXP);
    Rcpp::traits::input_parameter< size_t >::type nw(nwSEXP);
    rcpp_result_gen = Rcpp::wrap(C_glcm_textures_helper(x, w2, n_levels, shift, metrics, na_rm, impute_corr, ni, nw));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GLCMTextures_C_make_glcm_counts", (DL_FUNC) &_GLCMTextures_C_make_glcm_counts, 4},
    {"_GLCMTextures_C_make_glcm", (DL_FUNC) &_GLCMTextures_C_make_glcm, 4},
    {"_GLCMTextures_C_GLSV", (DL_FUNC) &_GLCMTextures_C_GLSV, 2},
    {"_GLCMTextures_C_glcm_metrics", (DL_FUNC) &_GLCMTextures_C_glcm_metrics, 7},
    {"_GLCMTextures_C_glcm_textures_helper", (DL_FUNC) &_GLCMTextures_C_glcm_textures_helper, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_GLCMTextures(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
