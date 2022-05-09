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
    Rcpp::IntegerMatrix X1C,
    Rcpp::IntegerMatrix a,
    Rcpp::IntegerMatrix d,    
    Rcpp::IntegerMatrix usg,
    Rcpp::IntegerMatrix usg_check,    
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
    s0 = X1C(a(k0, t0), t0) - 1; // ## this symbol to consider, 0-based
    int match_start = d(k0, t0);
    for(i = 0; i < St; i++) {
      if (match_start > pqs(i)) {
	pqs(i) = match_start;
      }
    }
    // now - where it goes - 0 based
    a(start_count[s0] + nso[s0], t0 + 1) = a(k0, t0);
    d(start_count[s0] + nso[s0], t0 + 1) = pqs(s0);
    usg.row(k0 + 1) = usg.row(k0);
    if (s0 < first_usg_minimal_symbol - 1) {
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
    for(s0 = 0; s0 < first_usg_minimal_symbol - 1; s0++) {
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






// from stackoverflow

// [[Rcpp::export]]
IntegerVector order_(IntegerVector x) {
  IntegerVector sorted = clone(x).sort();
  return match(sorted, x);
}



//' @export
// [[Rcpp::export]]
Rcpp::List Rcpp_ms_BuildIndices_Algorithm5(
    Rcpp::IntegerMatrix X1C,
    Rcpp::List all_symbols,
    bool verbose = false,
    bool do_checks = false,
    bool check_vs_indices = false,
    int egs = 100,
    int n_min_symboils = 100,
    bool with_Rcpp = true
) {    
    //
    const int K = X1C.nrow();
    const int T = X1C.ncol();
    //
    Rcpp::IntegerVector n_symbols_per_grid(all_symbols.length());
    int Smax = 0;
    int i, k;
    for(i = 0; i < all_symbols.length(); i++) {
      Rcpp::NumericMatrix temp = Rcpp::as<Rcpp::NumericMatrix>(all_symbols[i]);
      n_symbols_per_grid[i] = temp.nrow();
      if (n_symbols_per_grid[i] > Smax) {
	Smax = n_symbols_per_grid[i];
      }
    }
    // build arrays including a, d and now u, v, c
    Rcpp::IntegerMatrix a(K, T + 1);
    Rcpp::IntegerMatrix c(T + 1);
    Rcpp::IntegerVector b(K);
    for(k = 0; k < K; k++) {
      a[k, 0] = k;
    }
    a(_, 1) = order_(X1C(_, 0)) - 1;
    //
    // d
    //
    Rcpp::IntegerMatrix d(K + 1, T + 1);
    for(int i = 0; i < T; i++) {
      d(0, i) = i + 1;
      d(K, i) = i + 1;
    }
    // argh, will try and get rid of
    Rcpp::IntegerVector d_vec(K + 1);
    d_vec(0) = 1;
    d_vec(K) = 1;
    for(k = 0; k < K + 1; k++) {
      d(k, 1) = 0;
    }
    for(k = 0; k <= (K - 1 - 1); k++) {
      if (X1C(a(k, 1), 0) != X1C(a(k + 1, 1), 0)) {
	d_vec(k + 1) = 1;
      }
    }
    d_vec(0) = 2;
    d_vec(K) = 2;
    d(_, 1) = d_vec;
    //
    // ARGH THIS IS TORTURE - remove dependency on d_vec then continue from here
    //
    Rcpp::List to_return = Rcpp::List::create(
					      Rcpp::Named("a") = a,
					      Rcpp::Named("n_symbols_per_grid") = n_symbols_per_grid
    );
    return(to_return);
}

