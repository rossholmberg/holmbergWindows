// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// idDub
NumericVector idDub(int i, NumericVector inputlat, NumericVector inputlon, NumericMatrix inputdata, NumericVector outputlat, NumericVector outputlon);
RcppExport SEXP holmbergWindows_idDub(SEXP iSEXP, SEXP inputlatSEXP, SEXP inputlonSEXP, SEXP inputdataSEXP, SEXP outputlatSEXP, SEXP outputlonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type inputlat(inputlatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type inputlon(inputlonSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type inputdata(inputdataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type outputlat(outputlatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type outputlon(outputlonSEXP);
    rcpp_result_gen = Rcpp::wrap(idDub(i, inputlat, inputlon, inputdata, outputlat, outputlon));
    return rcpp_result_gen;
END_RCPP
}
