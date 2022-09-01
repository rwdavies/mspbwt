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
int Rcpp_decode_maximal_value_of_u(Rcpp::NumericMatrix& out_mat, Rcpp::NumericVector& out_vec, int v, int egs, bool do_checks);
RcppExport SEXP _mspbwt_Rcpp_decode_maximal_value_of_u(SEXP out_matSEXP, SEXP out_vecSEXP, SEXP vSEXP, SEXP egsSEXP, SEXP do_checksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type out_mat(out_matSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type out_vec(out_vecSEXP);
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type egs(egsSEXP);
    Rcpp::traits::input_parameter< bool >::type do_checks(do_checksSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_decode_maximal_value_of_u(out_mat, out_vec, v, egs, do_checks));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_encode_minimal_column_of_u
Rcpp::IntegerVector Rcpp_encode_minimal_column_of_u(Rcpp::IntegerVector& u);
RcppExport SEXP _mspbwt_Rcpp_encode_minimal_column_of_u(SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_encode_minimal_column_of_u(u));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_decode_minimal_value_of_u
int Rcpp_decode_minimal_value_of_u(Rcpp::IntegerVector& x, int v);
RcppExport SEXP _mspbwt_Rcpp_decode_minimal_value_of_u(SEXP xSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_decode_minimal_value_of_u(x, v));
    return rcpp_result_gen;
END_RCPP
}
// is_list
int is_list(SEXP x);
RcppExport SEXP _mspbwt_is_list(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(is_list(x));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_decode_value_of_usge
int Rcpp_decode_value_of_usge(Rcpp::List usge, int s, int v, int egs, int n_min_symbols);
RcppExport SEXP _mspbwt_Rcpp_decode_value_of_usge(SEXP usgeSEXP, SEXP sSEXP, SEXP vSEXP, SEXP egsSEXP, SEXP n_min_symbolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type usge(usgeSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type egs(egsSEXP);
    Rcpp::traits::input_parameter< int >::type n_min_symbols(n_min_symbolsSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_decode_value_of_usge(usge, s, v, egs, n_min_symbols));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_one_move_forward_buildindices
Rcpp::List Rcpp_one_move_forward_buildindices(Rcpp::IntegerMatrix X1C, Rcpp::IntegerMatrix a, Rcpp::IntegerMatrix d, Rcpp::IntegerMatrix usg, Rcpp::IntegerMatrix usg_check, int t, int K, Rcpp::IntegerVector symbol_count, int egs, int St, int n_min_symbols, bool do_checks);
RcppExport SEXP _mspbwt_Rcpp_one_move_forward_buildindices(SEXP X1CSEXP, SEXP aSEXP, SEXP dSEXP, SEXP usgSEXP, SEXP usg_checkSEXP, SEXP tSEXP, SEXP KSEXP, SEXP symbol_countSEXP, SEXP egsSEXP, SEXP StSEXP, SEXP n_min_symbolsSEXP, SEXP do_checksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type X1C(X1CSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type a(aSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type d(dSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type usg(usgSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type usg_check(usg_checkSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type symbol_count(symbol_countSEXP);
    Rcpp::traits::input_parameter< int >::type egs(egsSEXP);
    Rcpp::traits::input_parameter< int >::type St(StSEXP);
    Rcpp::traits::input_parameter< int >::type n_min_symbols(n_min_symbolsSEXP);
    Rcpp::traits::input_parameter< bool >::type do_checks(do_checksSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_one_move_forward_buildindices(X1C, a, d, usg, usg_check, t, K, symbol_count, egs, St, n_min_symbols, do_checks));
    return rcpp_result_gen;
END_RCPP
}
// order_
IntegerVector order_(IntegerVector x);
RcppExport SEXP _mspbwt_order_(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(order_(x));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_ms_BuildIndices_Algorithm5
Rcpp::List Rcpp_ms_BuildIndices_Algorithm5(Rcpp::IntegerMatrix X1C, Rcpp::List all_symbols, Rcpp::List indices, bool verbose, bool do_checks, bool check_vs_indices, int egs, int n_min_symbols, bool with_Rcpp);
RcppExport SEXP _mspbwt_Rcpp_ms_BuildIndices_Algorithm5(SEXP X1CSEXP, SEXP all_symbolsSEXP, SEXP indicesSEXP, SEXP verboseSEXP, SEXP do_checksSEXP, SEXP check_vs_indicesSEXP, SEXP egsSEXP, SEXP n_min_symbolsSEXP, SEXP with_RcppSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type X1C(X1CSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type all_symbols(all_symbolsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type do_checks(do_checksSEXP);
    Rcpp::traits::input_parameter< bool >::type check_vs_indices(check_vs_indicesSEXP);
    Rcpp::traits::input_parameter< int >::type egs(egsSEXP);
    Rcpp::traits::input_parameter< int >::type n_min_symbols(n_min_symbolsSEXP);
    Rcpp::traits::input_parameter< bool >::type with_Rcpp(with_RcppSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_ms_BuildIndices_Algorithm5(X1C, all_symbols, indices, verbose, do_checks, check_vs_indices, egs, n_min_symbols, with_Rcpp));
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
    {"_mspbwt_is_list", (DL_FUNC) &_mspbwt_is_list, 1},
    {"_mspbwt_Rcpp_decode_value_of_usge", (DL_FUNC) &_mspbwt_Rcpp_decode_value_of_usge, 5},
    {"_mspbwt_Rcpp_one_move_forward_buildindices", (DL_FUNC) &_mspbwt_Rcpp_one_move_forward_buildindices, 12},
    {"_mspbwt_order_", (DL_FUNC) &_mspbwt_order_, 1},
    {"_mspbwt_Rcpp_ms_BuildIndices_Algorithm5", (DL_FUNC) &_mspbwt_Rcpp_ms_BuildIndices_Algorithm5, 9},
    {"_mspbwt_BuildIndices_Algorithm5_Rcpp", (DL_FUNC) &_mspbwt_BuildIndices_Algorithm5_Rcpp, 1},
    {"_mspbwt_MatchZ_Algorithm5_Rcpp", (DL_FUNC) &_mspbwt_MatchZ_Algorithm5_Rcpp, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_mspbwt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
