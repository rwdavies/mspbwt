# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @export
Rcpp_encode_maximal_column_of_u <- function(u, egs, efficient = TRUE) {
    .Call('_mspbwt_Rcpp_encode_maximal_column_of_u', PACKAGE = 'mspbwt', u, egs, efficient)
}

#' @export
Rcpp_decode_maximal_value_of_u <- function(out_mat, out_vec, v, egs, do_checks = FALSE) {
    .Call('_mspbwt_Rcpp_decode_maximal_value_of_u', PACKAGE = 'mspbwt', out_mat, out_vec, v, egs, do_checks)
}

#' @export
Rcpp_encode_minimal_column_of_u <- function(u) {
    .Call('_mspbwt_Rcpp_encode_minimal_column_of_u', PACKAGE = 'mspbwt', u)
}

#' @export
Rcpp_decode_minimal_value_of_u <- function(x, v) {
    .Call('_mspbwt_Rcpp_decode_minimal_value_of_u', PACKAGE = 'mspbwt', x, v)
}

is_list <- function(x) {
    .Call('_mspbwt_is_list', PACKAGE = 'mspbwt', x)
}

#' @export
Rcpp_decode_value_of_usge <- function(usge, s, v, egs) {
    .Call('_mspbwt_Rcpp_decode_value_of_usge', PACKAGE = 'mspbwt', usge, s, v, egs)
}

#' @export
Rcpp_decode_value_of_usge_v2 <- function(usge, s, v, egs) {
    .Call('_mspbwt_Rcpp_decode_value_of_usge_v2', PACKAGE = 'mspbwt', usge, s, v, egs)
}

#' @export
Rcpp_get_k_given_encoded_u <- function(s, v2, usge, C, K, egs, verbose = FALSE) {
    .Call('_mspbwt_Rcpp_get_k_given_encoded_u', PACKAGE = 'mspbwt', s, v2, usge, C, K, egs, verbose)
}

#' @export
rcpp_go_backwards_one_step <- function(g, v, C, usge, egs, K, verbose = FALSE) {
    .Call('_mspbwt_rcpp_go_backwards_one_step', PACKAGE = 'mspbwt', g, v, C, usge, egs, K, verbose)
}

#' @export
Rcpp_find_index_backward <- function(g_in, v_in, all_symbols, usge_all, egs, K, list_of_columns_of_A, verbose = FALSE, use_list_of_columns_of_A = FALSE) {
    .Call('_mspbwt_Rcpp_find_index_backward', PACKAGE = 'mspbwt', g_in, v_in, all_symbols, usge_all, egs, K, list_of_columns_of_A, verbose, use_list_of_columns_of_A)
}

#' @export
Rcpp_get_f_given_Z <- function(Z, all_symbols, usge_all, egs, verbose = FALSE) {
    .Call('_mspbwt_Rcpp_get_f_given_Z', PACKAGE = 'mspbwt', Z, all_symbols, usge_all, egs, verbose)
}

#' @export
Rcpp_find_good_matches_without_a <- function(Z, all_symbols, usge_all, egs, pbwtL, pbwtM, K, hapMatcherR, which_snps_in_hapMatcherR, list_of_columns_of_A, use_list_of_columns_of_A, verbose = FALSE) {
    .Call('_mspbwt_Rcpp_find_good_matches_without_a', PACKAGE = 'mspbwt', Z, all_symbols, usge_all, egs, pbwtL, pbwtM, K, hapMatcherR, which_snps_in_hapMatcherR, list_of_columns_of_A, use_list_of_columns_of_A, verbose)
}

#' @export
Rcpp_one_move_forward_buildindices <- function(X1C, a, d, usg, usg_check, t, K, symbol_count, egs, St, n_min_symbols, do_checks) {
    .Call('_mspbwt_Rcpp_one_move_forward_buildindices', PACKAGE = 'mspbwt', X1C, a, d, usg, usg_check, t, K, symbol_count, egs, St, n_min_symbols, do_checks)
}

order_ <- function(x) {
    .Call('_mspbwt_order_', PACKAGE = 'mspbwt', x)
}

#' @export
Rcpp_ms_BuildIndices_Algorithm5 <- function(X1C, all_symbols, indices, verbose = FALSE, do_checks = FALSE, check_vs_indices = FALSE, egs = 100L, n_min_symbols = 100L, with_Rcpp = TRUE) {
    .Call('_mspbwt_Rcpp_ms_BuildIndices_Algorithm5', PACKAGE = 'mspbwt', X1C, all_symbols, indices, verbose, do_checks, check_vs_indices, egs, n_min_symbols, with_Rcpp)
}

#' @export
rcpp_wf <- function(k, t, s, usge_all, all_symbols, egs) {
    .Call('_mspbwt_rcpp_wf', PACKAGE = 'mspbwt', k, t, s, usge_all, all_symbols, egs)
}

#' @export
rcpp_find_restart <- function(e1, f1, g1, ec, fc, gc, X, XR, use_XR, a, Z, t, d, cols_to_use0, top_matches_list, K, verbose, use_cols_to_use0, test_d, use_d, cap_scan_count) {
    invisible(.Call('_mspbwt_rcpp_find_restart', PACKAGE = 'mspbwt', e1, f1, g1, ec, fc, gc, X, XR, use_XR, a, Z, t, d, cols_to_use0, top_matches_list, K, verbose, use_cols_to_use0, test_d, use_d, cap_scan_count))
}

#' @export
Rcpp_ms_MatchZ_Algorithm5 <- function(X, XR, ms_indices, Z, cols_to_use0, use_XR = FALSE, verbose = FALSE, do_checks = FALSE, check_vs_indices = FALSE, indices = FALSE, use_cols_to_use0 = FALSE, min_length = -1L, do_up_and_down_scan = FALSE, mspbwtL = 3L, mspbwtM = 3L, test_d = FALSE, have_d = TRUE, cap_scan_count = -1L) {
    .Call('_mspbwt_Rcpp_ms_MatchZ_Algorithm5', PACKAGE = 'mspbwt', X, XR, ms_indices, Z, cols_to_use0, use_XR, verbose, do_checks, check_vs_indices, indices, use_cols_to_use0, min_length, do_up_and_down_scan, mspbwtL, mspbwtM, test_d, have_d, cap_scan_count)
}

#' @export
BuildIndices_Algorithm5_Rcpp <- function(X) {
    .Call('_mspbwt_BuildIndices_Algorithm5_Rcpp', PACKAGE = 'mspbwt', X)
}

#' @export
MatchZ_Algorithm5_Rcpp <- function(X, a, u, v, c, d, Z, verbose = FALSE, do_checks = FALSE, start_rows = 100L) {
    .Call('_mspbwt_MatchZ_Algorithm5_Rcpp', PACKAGE = 'mspbwt', X, a, u, v, c, d, Z, verbose, do_checks, start_rows)
}

