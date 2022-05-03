#include <Rcpp.h>
using namespace Rcpp;



Rcpp::List Rcpp_encode_maximal_column_of_u(
    Rcpp::IntegerVector u,
    int egs,
    bool efficient = true
);


// likely to change this function substantially!

//' @export
// [[Rcpp::export]]
Rcpp::List Rcpp_one_move_forward_buildindices(
    Rcpp::NumericMatrix X1C,
    Rcpp::IntegerMatrix a,
    Rcpp::IntegerMatrix usg,
    Rcpp::IntegerMatrix usg_check,    
    Rcpp::IntegerVector d_vec,
    Rcpp::IntegerVector prev_d,
    int t,
    int K,
    Rcpp::IntegerVector symbol_count,
    int egs,
    int St,
    int n_min_symbols,
    bool do_checks
) {
  //
  int t0 = t - 1;
  int i;
  //
  Rcpp::List usge(St);
  int first_usg_minimal_symbol = 1; // 1-based
  int s0;
  for(s0 = 0; s0 < St; s0++) {
    if (symbol_count(s0) > n_min_symbols) {
      first_usg_minimal_symbol++;
    } else {
      Rcpp::IntegerVector temp_vec(symbol_count(s0));
      temp_vec.fill(-1);
      usge[s0] = temp_vec;
    }
  }
  //
  //
  //
  Rcpp::IntegerVector start_count(symbol_count.length() + 1);
  start_count(0) = 0;
  for(i = 0; i < symbol_count.length(); i++) {
    start_count(i + 1) = start_count(i) + symbol_count(i);
  }
  //
  //
  //
  Rcpp::IntegerVector nso(St);
  nso.fill(0);
  Rcpp::IntegerVector pqs(St);
  pqs.fill(t);
  usg.fill(0);
  //
  //
  //
  for(int k0 = 0; k0 < K; k0++) { // haps (1-based) (??)
    if (a(k0, t0) < 0) {
    std::cout << "k0 = " << k0 << std::endl;;
    std::cout << "t0 = " << t0 << std::endl;;    
    std::cout << "a(k0, t0) = " << a(k0, t0) << std::endl; // am here, have weird 0 result
    }
    s0 = X1C(a(k0, t0), t0) - 1; // ## this symbol to consider, 0-based
    int match_start = prev_d(k0);
    for(i = 0; i < St; i++) {
      if (match_start > pqs(i)) {
	pqs(i) = match_start;
      }
    }
    // now - where it goes - 0 based
    a(start_count[s0] + nso[s0], t0 + 1) = a(k0, t0);
    d_vec(start_count[s0] + nso[s0]) = pqs(s0);
    usg.row(k0 + 1) = usg.row(k0);
    if (s0 < first_usg_minimal_symbol) {
      usg(k0 + 1, s0) = usg(k0 + 1, s0) + 1;
    } else {
      Rcpp::IntegerVector temp = usge[s0];
      temp[nso[s0]] = k0 + 1;
      usge[s0] = temp;
    }
    pqs(s0) = 0;
    nso(s0)++;
    if (do_checks) {
      usg_check.row(k0 + 1) = usg_check.row(k0);
      usg_check(k0 + 1, s0)++;
    }
  }
  if ((first_usg_minimal_symbol - 1 - 1) >= 0) {
    for(s0 = 0; s0 < first_usg_minimal_symbol - 1 - 1; s0++) {
      usge[s0] = Rcpp_encode_maximal_column_of_u(usg(_, s0), egs);
    }
  }
  // std::cout << "d_vec = ";
  // for(int q = 100; q < 110; q++) {
  //   std::cout << d_vec(q) << ", ";
  // }
  Rcpp::List to_return = Rcpp::List::create(
	Rcpp::Named("usge") = usge
    );
  return(to_return);
}
