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
  int prev_value = symbol_count(0);
  for(s0 = 0; s0 < St; s0++) {
    if ((symbol_count(s0) > n_min_symbols) & (symbol_count(s0) <= prev_value)) {
      first_usg_minimal_symbol++;
    } else {
      Rcpp::IntegerVector temp_vec(symbol_count(s0));
      temp_vec.fill(-1);
      usge[s0] = temp_vec;
    }
    prev_value = symbol_count(s0);
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
    if (s0 == (-1)) {
      s0 = St - 1;
    }
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
      // this the problem with s0
      // what triggers this
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
    Rcpp::List indices,
    bool verbose = false,
    bool do_checks = false,
    bool check_vs_indices = false,
    int egs = 100,
    int n_min_symbols = 100,
    bool with_Rcpp = true
) {    
    //
    if (verbose) {
        std::cout << "initialize" << std::endl;
    }
    const int K = X1C.nrow();
    const int T = X1C.ncol();
    //
    // don't use this, but need to pass through
    //
    Rcpp::IntegerMatrix usg_check(1, 1);
    //
    //
    Rcpp::IntegerVector n_symbols_per_grid(all_symbols.length());
    int Smax = 0;
    int i, k, t;
    for(i = 0; i < all_symbols.length(); i++) {
      Rcpp::NumericMatrix temp = Rcpp::as<Rcpp::NumericMatrix>(all_symbols[i]);
      n_symbols_per_grid[i] = temp.nrow();
      if (n_symbols_per_grid[i] > Smax) {
	Smax = n_symbols_per_grid[i];
      }
    }
    //
    // extract this to change 0 to max val
    //
    Rcpp::IntegerVector X1C_col = X1C(_, 0);
    for(k = 0; k < K; k++) {
      if (X1C_col(k) == 0) {
	X1C_col(k)=n_symbols_per_grid(0);
      }
    }
    // build arrays including a, d and now u, v, c
    Rcpp::IntegerMatrix a(K, T + 1);
    Rcpp::IntegerMatrix c(T + 1);
    Rcpp::IntegerVector b(K);
    for(k = 0; k < K; k++) {
      a(k, 0) = k;
    }
    a(_, 1) = order_(X1C_col) - 1;
    //
    // d
    //
    Rcpp::IntegerMatrix d(K + 1, T + 1);
    d.fill(0); // unnecessary?
    for(k = 0; k <= (K - 1 - 1); k++) {
      if (X1C_col(a(k, 1)) != X1C_col(a(k + 1, 1))) {
	d(k + 1, 1) = 1;
      }
    }
    for(i = 0; i <= T; i++) {
      d(0, i) = i + 1;
      d(K, i) = i + 1;
    }
    //
    //
    //
    Rcpp::IntegerVector ns_obs(Smax);
    //
    Rcpp::List usge_all(T);
    int Smax_for_usl = 0;
    for(t = 0; t < T; t++) {
        Rcpp::NumericMatrix temp2 = Rcpp::as<Rcpp::NumericMatrix>(all_symbols[t]);
	int x = 0;
	for(i = 0; i < temp2.nrow(); i++) {
	    if (temp2(i, 1) > n_min_symbols) {
	        x += 1;
	    }
	}
	if (x > Smax_for_usl) {
	    Smax_for_usl = x;
	}
    }
    Rcpp::IntegerMatrix usg(K + 1, Smax_for_usl);
    //
    //
    //
    for(t = 1; t <= T; t++) {
      //
      //
      int St = n_symbols_per_grid(t - 1);
      Rcpp::NumericMatrix temp2 = Rcpp::as<Rcpp::NumericMatrix>(all_symbols[t - 1]);
      // argh
      Rcpp::IntegerVector symbol_count(temp2.nrow());
      for(i = 0; i < temp2.nrow(); i++) {
	symbol_count(i) = temp2(i, 1);
      }
      Rcpp::List usge = Rcpp_one_move_forward_buildindices(
          X1C, a, d, usg, usg_check, t, K, symbol_count, egs, St, n_min_symbols, do_checks
      );
      usge_all(t - 1) = usge("usge");
      if (t == 1) {
	d(0, t) = t + 1;
      }
    }
    Rcpp::List to_return = Rcpp::List::create(
        Rcpp::Named("a") = a,
        Rcpp::Named("d") = d,
        Rcpp::Named("usge_all") = usge_all,
	Rcpp::Named("egs") = egs,
	Rcpp::Named("n_min_symbols") = n_min_symbols,
	Rcpp::Named("all_symbols") = all_symbols
    );
    return(to_return);
}







// rcpp_wf(k, t, s, usge_all, all_symbols, egs) {
//     if (s == 0) {
//         Rcpp::NumericMatrix temp_mat = all_symbols[t];
//         s = temp_mat.nrow();
//     }
//     int u = Rcpp_decode_value_of_usge(usge_all[t - 1], s, k, egs, n_min_symbols);
    
	  
//     c <- c(0, cumsum(all_symbols[[t]][, 2]))[s]

//     ) + c
//     u
// }


  

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix Rcpp_ms_MatchZ_Algorithm5(
    Rcpp::NumericMatrix X,
    Rcpp::List ms_indices,
    Rcpp::IntegerVector Z,
    bool verbose = false,
    bool do_checks  = false,
    bool check_vs_indices = false,
    bool indices = false
) {
    int K = X.nrow();
    int T = X.ncol();
    int k;
    //
    // get things out of lists so we can use them
    //
    Rcpp::NumericMatrix a = ms_indices["a"];
    Rcpp::NumericMatrix d = ms_indices["d"];
    Rcpp::List usge_all = ms_indices["usge_all"];
    int egs = ms_indices["egs"];
    int n_min_symbols = ms_indices["n_min_symbols"];
    Rcpp::List all_symbols = ms_indices["all_symbols"];
    //
    // initialize
    //
    Rcpp::IntegerVector e(T);
    Rcpp::IntegerVector f(T);
    Rcpp::IntegerVector g(T);    
    e(0) = 0;

    Rcpp::NumericMatrix temp_mat = all_symbols[0];
    Rcpp::NumericVector x(temp_mat.nrow() + 1);
    x(0) = 0;
    for(int i = 1; i < x.length(); i++) {
      x(i) = x(i) + temp_mat(i - 1, 1);
    }
    int Z_1 = Z(0);
    if (Z_1 == 0) {
      Z_1 = temp_mat.nrow();
    }
    f(0) = x(Z_1 - 1);
    g(0) = x(Z_1);
    //
    // just do easy bit for now
    //
    int fc = f(0);
    int gc = g(0);
    int ec = e(0);
    int e1 = -1;
    Rcpp::List top_matches_list; // probably fine unless this becomes massive!
    int i_top_matches = 0;
    //
    // loop city
    //
    // t stays 1-BASED
    for(int t = 2; t < T; t++) {
      int f1, g1;
      //f1 = wf(fc, t, Z[t - 1], usge_all, all_symbols, egs);
      //g1 = wf(gc, t, Z[t - 1], usge_all, all_symbols, egs);
    }
    Rcpp::NumericMatrix top_matches(5, 4);
    return(top_matches);
}
