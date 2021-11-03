// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Rcpp_encode_maximal_column_of_u
Rcpp::List Rcpp_encode_maximal_column_of_u(Rcpp::IntegerVector u, int egs, bool efficient);
RcppExport SEXP _mspbwt_Rcpp_encode_maximal_column_of_u(SEXP uSEXP, SEXP egsSEXP, SEXP efficientSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< int >::type egs(egsSEXP);
    Rcpp::traits::input_parameter< bool >::type efficient(efficientSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_encode_maximal_column_of_u(u, egs, efficient));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_decode_maximal_value_of_u
int Rcpp_decode_maximal_value_of_u(Rcpp::NumericMatrix out_mat, Rcpp::NumericVector out_vec, int v, int egs, bool do_checks);
RcppExport SEXP _mspbwt_Rcpp_decode_maximal_value_of_u(SEXP out_matSEXP, SEXP out_vecSEXP, SEXP vSEXP, SEXP egsSEXP, SEXP do_checksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type out_mat(out_matSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type out_vec(out_vecSEXP);
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type egs(egsSEXP);
    Rcpp::traits::input_parameter< bool >::type do_checks(do_checksSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_decode_maximal_value_of_u(out_mat, out_vec, v, egs, do_checks));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_encode_minimal_column_of_u
Rcpp::IntegerVector Rcpp_encode_minimal_column_of_u(Rcpp::IntegerVector u);
RcppExport SEXP _mspbwt_Rcpp_encode_minimal_column_of_u(SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_encode_minimal_column_of_u(u));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_decode_minimal_value_of_u
int Rcpp_decode_minimal_value_of_u(Rcpp::IntegerVector x, int v);
RcppExport SEXP _mspbwt_Rcpp_decode_minimal_value_of_u(SEXP xSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_decode_minimal_value_of_u(x, v));
    return rcpp_result_gen;
END_RCPP
}
// BuildIndices_Algorithm5_Rcpp
Rcpp::List BuildIndices_Algorithm5_Rcpp(Rcpp::IntegerMatrix& X);
RcppExport SEXP _mspbwt_BuildIndices_Algorithm5_Rcpp(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(BuildIndices_Algorithm5_Rcpp(X));
    return rcpp_result_gen;
END_RCPP
}
// MatchZ_Algorithm5_Rcpp
Rcpp::NumericMatrix MatchZ_Algorithm5_Rcpp(Rcpp::IntegerMatrix& X, Rcpp::IntegerMatrix& a, Rcpp::NumericMatrix& u, Rcpp::NumericMatrix& v, Rcpp::NumericVector& c, Rcpp::NumericMatrix& d, Rcpp::NumericVector& Z, bool verbose, bool do_checks, int start_rows);
RcppExport SEXP _mspbwt_MatchZ_Algorithm5_Rcpp(SEXP XSEXP, SEXP aSEXP, SEXP uSEXP, SEXP vSEXP, SEXP cSEXP, SEXP dSEXP, SEXP ZSEXP, SEXP verboseSEXP, SEXP do_checksSEXP, SEXP start_rowsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix& >::type a(aSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type u(uSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type v(vSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type c(cSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type d(dSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type do_checks(do_checksSEXP);
    Rcpp::traits::input_parameter< int >::type start_rows(start_rowsSEXP);
    rcpp_result_gen = Rcpp::wrap(MatchZ_Algorithm5_Rcpp(X, a, u, v, c, d, Z, verbose, do_checks, start_rows));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mspbwt_Rcpp_encode_maximal_column_of_u", (DL_FUNC) &_mspbwt_Rcpp_encode_maximal_column_of_u, 3},
    {"_mspbwt_Rcpp_decode_maximal_value_of_u", (DL_FUNC) &_mspbwt_Rcpp_decode_maximal_value_of_u, 5},
    {"_mspbwt_Rcpp_encode_minimal_column_of_u", (DL_FUNC) &_mspbwt_Rcpp_encode_minimal_column_of_u, 1},
    {"_mspbwt_Rcpp_decode_minimal_value_of_u", (DL_FUNC) &_mspbwt_Rcpp_decode_minimal_value_of_u, 2},
    {"_mspbwt_BuildIndices_Algorithm5_Rcpp", (DL_FUNC) &_mspbwt_BuildIndices_Algorithm5_Rcpp, 1},
    {"_mspbwt_MatchZ_Algorithm5_Rcpp", (DL_FUNC) &_mspbwt_MatchZ_Algorithm5_Rcpp, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_mspbwt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
