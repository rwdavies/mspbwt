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
int Rcpp_decode_maximal_value_of_u(Rcpp::IntegerMatrix& out_mat, Rcpp::IntegerVector& out_vec, int v, int egs, bool do_checks);
RcppExport SEXP _mspbwt_Rcpp_decode_maximal_value_of_u(SEXP out_matSEXP, SEXP out_vecSEXP, SEXP vSEXP, SEXP egsSEXP, SEXP do_checksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix& >::type out_mat(out_matSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type out_vec(out_vecSEXP);
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
int Rcpp_decode_value_of_usge(Rcpp::List usge, int s, int v, int egs);
RcppExport SEXP _mspbwt_Rcpp_decode_value_of_usge(SEXP usgeSEXP, SEXP sSEXP, SEXP vSEXP, SEXP egsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type usge(usgeSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type egs(egsSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_decode_value_of_usge(usge, s, v, egs));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_decode_value_of_usge_v2
int Rcpp_decode_value_of_usge_v2(Rcpp::List& usge, int s, int v, int egs);
RcppExport SEXP _mspbwt_Rcpp_decode_value_of_usge_v2(SEXP usgeSEXP, SEXP sSEXP, SEXP vSEXP, SEXP egsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type usge(usgeSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type egs(egsSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_decode_value_of_usge_v2(usge, s, v, egs));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_get_k_given_encoded_u
int Rcpp_get_k_given_encoded_u(int s, int v2, Rcpp::List& usge, Rcpp::IntegerMatrix& C, int K, int egs, bool verbose);
RcppExport SEXP _mspbwt_Rcpp_get_k_given_encoded_u(SEXP sSEXP, SEXP v2SEXP, SEXP usgeSEXP, SEXP CSEXP, SEXP KSEXP, SEXP egsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type v2(v2SEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type usge(usgeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix& >::type C(CSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type egs(egsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_get_k_given_encoded_u(s, v2, usge, C, K, egs, verbose));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_go_backwards_one_step
int rcpp_go_backwards_one_step(int g, int v, Rcpp::IntegerMatrix& C, Rcpp::List& usge, int egs, int K, bool verbose);
RcppExport SEXP _mspbwt_rcpp_go_backwards_one_step(SEXP gSEXP, SEXP vSEXP, SEXP CSEXP, SEXP usgeSEXP, SEXP egsSEXP, SEXP KSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type g(gSEXP);
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix& >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type usge(usgeSEXP);
    Rcpp::traits::input_parameter< int >::type egs(egsSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_go_backwards_one_step(g, v, C, usge, egs, K, verbose));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_find_index_backward
int Rcpp_find_index_backward(int g_in, int v_in, Rcpp::List& all_symbols, Rcpp::List& usge_all, int egs, int K, Rcpp::List& list_of_columns_of_A, bool verbose, bool use_list_of_columns_of_A);
RcppExport SEXP _mspbwt_Rcpp_find_index_backward(SEXP g_inSEXP, SEXP v_inSEXP, SEXP all_symbolsSEXP, SEXP usge_allSEXP, SEXP egsSEXP, SEXP KSEXP, SEXP list_of_columns_of_ASEXP, SEXP verboseSEXP, SEXP use_list_of_columns_of_ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type g_in(g_inSEXP);
    Rcpp::traits::input_parameter< int >::type v_in(v_inSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type all_symbols(all_symbolsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type usge_all(usge_allSEXP);
    Rcpp::traits::input_parameter< int >::type egs(egsSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type list_of_columns_of_A(list_of_columns_of_ASEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type use_list_of_columns_of_A(use_list_of_columns_of_ASEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_find_index_backward(g_in, v_in, all_symbols, usge_all, egs, K, list_of_columns_of_A, verbose, use_list_of_columns_of_A));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_find_good_matches_without_a
Rcpp::IntegerMatrix Rcpp_find_good_matches_without_a(Rcpp::IntegerVector& f, Rcpp::IntegerVector& Z, Rcpp::List& all_symbols, Rcpp::List& usge_all, int egs, int pbwtL, int pbwtM, int K, Rcpp::RawMatrix& hapMatcherR, Rcpp::IntegerVector& which_snps_in_hapMatcherR, Rcpp::List& list_of_columns_of_A, bool use_list_of_columns_of_A, bool verbose);
RcppExport SEXP _mspbwt_Rcpp_find_good_matches_without_a(SEXP fSEXP, SEXP ZSEXP, SEXP all_symbolsSEXP, SEXP usge_allSEXP, SEXP egsSEXP, SEXP pbwtLSEXP, SEXP pbwtMSEXP, SEXP KSEXP, SEXP hapMatcherRSEXP, SEXP which_snps_in_hapMatcherRSEXP, SEXP list_of_columns_of_ASEXP, SEXP use_list_of_columns_of_ASEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type f(fSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type all_symbols(all_symbolsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type usge_all(usge_allSEXP);
    Rcpp::traits::input_parameter< int >::type egs(egsSEXP);
    Rcpp::traits::input_parameter< int >::type pbwtL(pbwtLSEXP);
    Rcpp::traits::input_parameter< int >::type pbwtM(pbwtMSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< Rcpp::RawMatrix& >::type hapMatcherR(hapMatcherRSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type which_snps_in_hapMatcherR(which_snps_in_hapMatcherRSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type list_of_columns_of_A(list_of_columns_of_ASEXP);
    Rcpp::traits::input_parameter< bool >::type use_list_of_columns_of_A(use_list_of_columns_of_ASEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_find_good_matches_without_a(f, Z, all_symbols, usge_all, egs, pbwtL, pbwtM, K, hapMatcherR, which_snps_in_hapMatcherR, list_of_columns_of_A, use_list_of_columns_of_A, verbose));
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
// rcpp_wf
int rcpp_wf(int k, int t, int s, Rcpp::List& usge_all, Rcpp::List& all_symbols, int egs);
RcppExport SEXP _mspbwt_rcpp_wf(SEXP kSEXP, SEXP tSEXP, SEXP sSEXP, SEXP usge_allSEXP, SEXP all_symbolsSEXP, SEXP egsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type usge_all(usge_allSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type all_symbols(all_symbolsSEXP);
    Rcpp::traits::input_parameter< int >::type egs(egsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_wf(k, t, s, usge_all, all_symbols, egs));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_find_restart
void rcpp_find_restart(int& e1, int& f1, int& g1, int& ec, int& fc, int& gc, Rcpp::IntegerMatrix& X, Rcpp::RawMatrix& XR, bool use_XR, Rcpp::IntegerMatrix& a, Rcpp::IntegerVector& Z, int& t, Rcpp::IntegerMatrix& d, Rcpp::IntegerVector& cols_to_use0, Rcpp::List& top_matches_list, int K, bool verbose, bool use_cols_to_use0, bool test_d, bool use_d, int cap_scan_count);
RcppExport SEXP _mspbwt_rcpp_find_restart(SEXP e1SEXP, SEXP f1SEXP, SEXP g1SEXP, SEXP ecSEXP, SEXP fcSEXP, SEXP gcSEXP, SEXP XSEXP, SEXP XRSEXP, SEXP use_XRSEXP, SEXP aSEXP, SEXP ZSEXP, SEXP tSEXP, SEXP dSEXP, SEXP cols_to_use0SEXP, SEXP top_matches_listSEXP, SEXP KSEXP, SEXP verboseSEXP, SEXP use_cols_to_use0SEXP, SEXP test_dSEXP, SEXP use_dSEXP, SEXP cap_scan_countSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int& >::type e1(e1SEXP);
    Rcpp::traits::input_parameter< int& >::type f1(f1SEXP);
    Rcpp::traits::input_parameter< int& >::type g1(g1SEXP);
    Rcpp::traits::input_parameter< int& >::type ec(ecSEXP);
    Rcpp::traits::input_parameter< int& >::type fc(fcSEXP);
    Rcpp::traits::input_parameter< int& >::type gc(gcSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::RawMatrix& >::type XR(XRSEXP);
    Rcpp::traits::input_parameter< bool >::type use_XR(use_XRSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix& >::type a(aSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int& >::type t(tSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix& >::type d(dSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type cols_to_use0(cols_to_use0SEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type top_matches_list(top_matches_listSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type use_cols_to_use0(use_cols_to_use0SEXP);
    Rcpp::traits::input_parameter< bool >::type test_d(test_dSEXP);
    Rcpp::traits::input_parameter< bool >::type use_d(use_dSEXP);
    Rcpp::traits::input_parameter< int >::type cap_scan_count(cap_scan_countSEXP);
    rcpp_find_restart(e1, f1, g1, ec, fc, gc, X, XR, use_XR, a, Z, t, d, cols_to_use0, top_matches_list, K, verbose, use_cols_to_use0, test_d, use_d, cap_scan_count);
    return R_NilValue;
END_RCPP
}
// Rcpp_ms_MatchZ_Algorithm5
Rcpp::NumericMatrix Rcpp_ms_MatchZ_Algorithm5(Rcpp::IntegerMatrix& X, Rcpp::RawMatrix& XR, Rcpp::List& ms_indices, Rcpp::IntegerVector& Z, Rcpp::IntegerVector& cols_to_use0, bool use_XR, bool verbose, bool do_checks, bool check_vs_indices, bool indices, bool use_cols_to_use0, int min_length, bool do_up_and_down_scan, int mspbwtL, int mspbwtM, bool test_d, bool have_d, int cap_scan_count);
RcppExport SEXP _mspbwt_Rcpp_ms_MatchZ_Algorithm5(SEXP XSEXP, SEXP XRSEXP, SEXP ms_indicesSEXP, SEXP ZSEXP, SEXP cols_to_use0SEXP, SEXP use_XRSEXP, SEXP verboseSEXP, SEXP do_checksSEXP, SEXP check_vs_indicesSEXP, SEXP indicesSEXP, SEXP use_cols_to_use0SEXP, SEXP min_lengthSEXP, SEXP do_up_and_down_scanSEXP, SEXP mspbwtLSEXP, SEXP mspbwtMSEXP, SEXP test_dSEXP, SEXP have_dSEXP, SEXP cap_scan_countSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::RawMatrix& >::type XR(XRSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type ms_indices(ms_indicesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type cols_to_use0(cols_to_use0SEXP);
    Rcpp::traits::input_parameter< bool >::type use_XR(use_XRSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type do_checks(do_checksSEXP);
    Rcpp::traits::input_parameter< bool >::type check_vs_indices(check_vs_indicesSEXP);
    Rcpp::traits::input_parameter< bool >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< bool >::type use_cols_to_use0(use_cols_to_use0SEXP);
    Rcpp::traits::input_parameter< int >::type min_length(min_lengthSEXP);
    Rcpp::traits::input_parameter< bool >::type do_up_and_down_scan(do_up_and_down_scanSEXP);
    Rcpp::traits::input_parameter< int >::type mspbwtL(mspbwtLSEXP);
    Rcpp::traits::input_parameter< int >::type mspbwtM(mspbwtMSEXP);
    Rcpp::traits::input_parameter< bool >::type test_d(test_dSEXP);
    Rcpp::traits::input_parameter< bool >::type have_d(have_dSEXP);
    Rcpp::traits::input_parameter< int >::type cap_scan_count(cap_scan_countSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_ms_MatchZ_Algorithm5(X, XR, ms_indices, Z, cols_to_use0, use_XR, verbose, do_checks, check_vs_indices, indices, use_cols_to_use0, min_length, do_up_and_down_scan, mspbwtL, mspbwtM, test_d, have_d, cap_scan_count));
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
    {"_mspbwt_Rcpp_decode_value_of_usge", (DL_FUNC) &_mspbwt_Rcpp_decode_value_of_usge, 4},
    {"_mspbwt_Rcpp_decode_value_of_usge_v2", (DL_FUNC) &_mspbwt_Rcpp_decode_value_of_usge_v2, 4},
    {"_mspbwt_Rcpp_get_k_given_encoded_u", (DL_FUNC) &_mspbwt_Rcpp_get_k_given_encoded_u, 7},
    {"_mspbwt_rcpp_go_backwards_one_step", (DL_FUNC) &_mspbwt_rcpp_go_backwards_one_step, 7},
    {"_mspbwt_Rcpp_find_index_backward", (DL_FUNC) &_mspbwt_Rcpp_find_index_backward, 9},
    {"_mspbwt_Rcpp_find_good_matches_without_a", (DL_FUNC) &_mspbwt_Rcpp_find_good_matches_without_a, 13},
    {"_mspbwt_Rcpp_one_move_forward_buildindices", (DL_FUNC) &_mspbwt_Rcpp_one_move_forward_buildindices, 12},
    {"_mspbwt_order_", (DL_FUNC) &_mspbwt_order_, 1},
    {"_mspbwt_Rcpp_ms_BuildIndices_Algorithm5", (DL_FUNC) &_mspbwt_Rcpp_ms_BuildIndices_Algorithm5, 9},
    {"_mspbwt_rcpp_wf", (DL_FUNC) &_mspbwt_rcpp_wf, 6},
    {"_mspbwt_rcpp_find_restart", (DL_FUNC) &_mspbwt_rcpp_find_restart, 21},
    {"_mspbwt_Rcpp_ms_MatchZ_Algorithm5", (DL_FUNC) &_mspbwt_Rcpp_ms_MatchZ_Algorithm5, 18},
    {"_mspbwt_BuildIndices_Algorithm5_Rcpp", (DL_FUNC) &_mspbwt_BuildIndices_Algorithm5_Rcpp, 1},
    {"_mspbwt_MatchZ_Algorithm5_Rcpp", (DL_FUNC) &_mspbwt_MatchZ_Algorithm5_Rcpp, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_mspbwt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
