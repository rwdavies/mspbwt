#include <Rcpp.h>
using namespace Rcpp;



Rcpp::List Rcpp_encode_maximal_column_of_u(
    Rcpp::IntegerVector u,
    int egs,
    bool efficient = true
);

int Rcpp_decode_value_of_usge(
    Rcpp::List usge,
    int s,
    int v,
    int egs
);

int Rcpp_decode_value_of_usge_v2(
    Rcpp::List & usge,
    int s,
    int v,
    int egs
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






//' @export
// [[Rcpp::export]]
int rcpp_wf(
    int k,
    int t,
    int s,
    Rcpp::List & usge_all,
    Rcpp::List & all_symbols,
    int egs
) {
    Rcpp::NumericMatrix temp_mat = Rcpp::as<Rcpp::NumericMatrix>(all_symbols[t - 1]);      
    // recall t is 1-based here
    if (s == 0) {
      s = temp_mat.nrow();
    }
    int u = Rcpp_decode_value_of_usge(usge_all[t - 1], s, k, egs);
    int c = 0;
    // if s is 1, we are good. everything else we add
    if (s > 1) {
      for(int s2 = 1; s2 < s; s2++) {
	c += temp_mat(s2 - 1, 1);
      }
    }
    u += c;
    return(u);
}














//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix Rcpp_ms_MatchZ_Algorithm5(
    Rcpp::IntegerMatrix & X,
    Rcpp::RawMatrix & XR,    
    Rcpp::List & ms_indices,
    Rcpp::IntegerVector & Z,
    Rcpp::IntegerVector & cols_to_use0,
    bool use_XR = false,
    bool verbose = false,
    bool do_checks  = false,
    bool check_vs_indices = false,
    bool indices = false,
    bool use_cols_to_use0 = false,
    int min_length = -1,
    bool do_up_and_down_scan = false,
    int mspbwtL = 3,
    int mspbwtM = 3
) {
    if (verbose) {
      std::cout << "Inside Rcpp ms algorithm" << std::endl;
    }
    int K;
    int T;
    if (use_XR) {
      K = XR.nrow();
      T = XR.ncol();            
    } else {
      K = X.nrow();
      T = X.ncol();      
    }
    if (use_cols_to_use0) {
      T = cols_to_use0.length();
    }    
    int k, index, t;
    int f1, g1, e1_local;
    bool matches_lower, matches_upper;
    //
    // get things out of lists so we can use them
    //
    Rcpp::IntegerMatrix a = ms_indices["a"];
    Rcpp::IntegerMatrix d = ms_indices["d"];
    Rcpp::List usge_all = ms_indices["usge_all"];
    int egs = ms_indices["egs"];
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
      x(i) = x(i - 1) + temp_mat(i - 1, 1);
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
    //
    // initialize up and down scan, potentially
    //
    // uppy downy = ud
    // a = above, b = below
    Rcpp::IntegerVector ud_up_prev(mspbwtL);
    Rcpp::IntegerVector ud_up_cur(mspbwtL);
    Rcpp::IntegerVector ud_up_length_prev(mspbwtL);
    Rcpp::IntegerVector ud_up_length_cur(mspbwtL);
    //
    Rcpp::IntegerVector ud_down_prev(mspbwtL);
    Rcpp::IntegerVector ud_down_cur(mspbwtL);
    Rcpp::IntegerVector ud_down_length_prev(mspbwtL);
    Rcpp::IntegerVector ud_down_length_cur(mspbwtL);
    int fg;
    int i0 = 1;
    int i0_cur, i0_prev, prev, cur, len;
    Rcpp::IntegerVector temp_vector(3);
    Rcpp::List uppy_downy_reporter;
    int uppy_downy_count = 0;
    // make
    Rcpp::IntegerMatrix uppy_downy_matrix(100, 4);    
    if (do_up_and_down_scan) {
        if (verbose) {
          std::cout << "Initialize up and down scan" << std::endl;
	}
        // 
        ud_up_prev.fill(-1);
	ud_down_prev.fill(-1);
	// initialize up
	// should automatically floor? as an integer
	fg = ((f(0) + g(0) - 1) / 2); // 0-based, include in "up"
	fg = f(0);
	// do up, include first entry fg
	i0 = 0; // ## 0-based
	while(i0 <= (mspbwtL - 1)) {
	  if (0 <= (fg - i0)) {
	    ud_up_prev(i0) = a(fg - i0, 1); // ## go up, so subtract
	    ud_up_length_prev(i0) = 0;
	  } else {
	    ud_up_prev(i0) = -1;
	    ud_up_length_prev(i0) = -1;
	  }
	  i0++;
	}
	// do down, go after first entry fg
	i0 = 0; // ## 0-based      
	while(i0 <= (mspbwtL - 1)) {
	  if ((fg + i0 + 1) <= (K - 1)) {
	    ud_down_prev(i0) = a(fg + i0 + 1, 1); // ## go up, so subtract
	    ud_down_length_prev(i0) = 0;
	  } else {
	    ud_down_prev(i0) = -1;
	    ud_down_length_prev(i0) = -1;
	  }
	  i0++;
	}
        if (verbose) {
          std::cout << "Done initialize up and down scan" << std::endl;
	}
    }
    i0 = -100; // should prevent problems?
    //
    // loop city
    //
    // t stays 1-BASED
    for(t = 2; t <= T; t++) {
      f1 = rcpp_wf(fc, t, Z(t - 1), usge_all, all_symbols, egs);
      g1 = rcpp_wf(gc, t, Z(t - 1), usge_all, all_symbols, egs);
      if (verbose) {
	std::cout << "Start of loop t=" <<  t << ", fc = " << fc << ",  gc = " << gc << ",  ec = " << ec << ",  Z[t - 1] = " << Z(t - 1) << ",  f1=" << f1 << ",  g1=" <<  g1 << ",  e1 = " << e1 << std::endl;
      }
      //
      // up and down scan bit
      //
      if (do_up_and_down_scan) {
	//fg = ((f1 + g1 - 1) / 2); // 0-based, include in "up"
	// not ideal but OK
	fg = f1; 
	//
	if (verbose) {
	  std::cout << "scan up" << std::endl;	// go "up" i.e. above i.e. up in the matrix
	}
	i0_cur =0; // 0-based, through local
	i0_prev = 0;
	// go through previous values
        while((i0_prev <= (mspbwtL - 1)) && (-1 < ud_up_prev(i0_cur)) && (0 <= (fg - i0_cur))) {
	    // focus on going through past list
	    prev = ud_up_prev(i0_prev);
	    // now what is the current, does that work
	    cur = a(fg - i0_cur, t); // ## go up, so subtract
	    if (cur == prev) {
	      // it is a match. save and increment match
	      ud_up_cur(i0_cur) = cur;
	      //## do the same for length. use previous length
	      ud_up_length_cur(i0_cur) = ud_up_length_prev(i0_prev) + 1;
	      //## increment up one
	      i0_cur++;
	    } else {
	      //## it is not a match, report it
	      //## do not increment cur
	      len = ud_up_length_prev(i0_prev);// ## 0-based
	      if (mspbwtM <= len) {
		// Rcpp::IntegerVector temp_vector(3);			  
		// temp_vector(0) = prev;
		// temp_vector(1) = t - 1;
		// temp_vector(2) = len;
		if (uppy_downy_count >= uppy_downy_matrix.nrow()) {
		  // double the size!
		  Rcpp::IntegerMatrix new_uppy_downy_matrix(uppy_downy_matrix.nrow() * 2, 3);
		  for(int i = 0; i < uppy_downy_matrix.nrow(); i++) {
		    new_uppy_downy_matrix(i, _) = uppy_downy_matrix(i, _);
		  }
		  uppy_downy_matrix=new_uppy_downy_matrix;
		}
		// if (prev == 542061 && len == (t - 1)) {
		//   std::cout << "Ahere I am, t = " << t << std::endl;
		//   std::cout << "uppy_downy_count = " << uppy_downy_count << std::endl;
		//   std::cout << "i0_cur = " << i0_cur << ", cur = " << cur << std::endl;
		//   std::cout << "i0_prev = " << i0_prev << ", prev = " << prev << std::endl;
		// }
		uppy_downy_matrix(uppy_downy_count, 0) = prev;
		uppy_downy_matrix(uppy_downy_count, 1) = t - 1;
		uppy_downy_matrix(uppy_downy_count, 2) = len;
		uppy_downy_count++;
		//uppy_downy_reporter.push_back(temp_vector);
	      }
	    }
	    i0_prev++;
	}
	if (verbose) {
	  std::cout << "fill in" << std::endl;	//## now fill in what was not set
	}
	while((i0_cur <= (mspbwtL - 1))) {
	    if (0 <= (fg - i0_cur)) {
	        cur = a(fg - i0_cur, t); // ## go up, so subtract
		ud_up_cur(i0_cur) = cur;
		// go backward, find start
		e1 = t; // recall t is 1-based
		if (!use_cols_to_use0) {
		  if (!use_XR) {
		      while((1 <= e1) && (X(cur, e1 - 1) == Z(e1 - 1))) {
			e1--;
		      }
		    } else {
		      while((1 <= e1) && (XR(cur, e1 - 1) == Z(e1 - 1))) {
			e1--;
		      }
		    }
		} else {
		  e1 = t;
		  if (!use_XR) {
		    while((1 <= e1) && (X(cur, cols_to_use0(e1 - 1)) == Z(e1 - 1))) {
		      e1--;
		    }
		  } else {
		    while((1 <= e1) && (XR(cur, cols_to_use0(e1 - 1)) == Z(e1 - 1))) {
		      e1--;
		    }
		  }
		}
		// go backwards, sort out start, how many before
		ud_up_length_cur(i0_cur) = t - e1;
	    } else {
	      ud_up_cur(i0_cur) = -1;
	      ud_up_length_cur(i0_cur) = -1;
	    }
	    i0_cur++;
	}
	// if (t <= 3) {
	//   std::cout << "----------------- t = " << t << " ------ RCPP" << std::endl;
	//   std::cout << "fg = " << fg << std::endl;
	//   std::cout << "ud_up_prev";
	//   std::cout << ud_up_prev;
	//   std::cout << std::endl;
	//   std::cout << "ud_up_length_prev";
	//   std::cout << ud_up_length_prev;
	//   std::cout << std::endl;
	//   std::cout << "ud_up_cur";
	//   std::cout << ud_up_cur;
	//   std::cout << std::endl;
	//   std::cout << "ud_up_length_cur";
	//   std::cout << ud_up_length_cur;
	//   std::cout << std::endl;
	// }
	ud_up_prev = ud_up_cur;
	ud_up_length_prev = ud_up_length_cur;
	if (verbose) {
	  std::cout << "scan down" << std::endl;	// go "down i.e. below in the matrix
	}
	i0_cur = 0; // ## 0-based, through local
	i0_prev = 0;
	//## go through previous values
	while((i0_prev <= (mspbwtL - 1)) && (-1 < ud_down_prev(i0_cur)) && ((fg + i0_cur + 1) <= (K - 1))) {
	  //## focus on going through past list
	  prev = ud_down_prev(i0_prev);
	  //## now what is the current, does that work
	  cur = a(fg + i0_cur + 1, t);
	  if (cur == prev) {
	    //## it is a match. save and increment match
	    ud_down_cur(i0_cur) = cur;
	    //## do the same for length. use previous length
	    ud_down_length_cur(i0_cur) = ud_down_length_prev(i0_prev) + 1;
	    //## increment up one
	    i0_cur++;
	  } else {
	    //## it is not a match, report it
	    //## do not increment cur
	    len = ud_down_length_prev(i0_prev); // ## 0-based
	    if (mspbwtM <= len) {
		// Rcpp::IntegerVector temp_vector(3);			  
		// temp_vector(0) = prev;
		// temp_vector(1) = t - 1;
		// temp_vector(2) = len;
		if (uppy_downy_count >= uppy_downy_matrix.nrow()) {
		  // double the size!
		  Rcpp::IntegerMatrix new_uppy_downy_matrix(uppy_downy_matrix.nrow() * 2, 3);
		  for(int i = 0; i < uppy_downy_matrix.nrow(); i++) {
		    new_uppy_downy_matrix(i, _) = uppy_downy_matrix(i, _);
		  }
		  uppy_downy_matrix=new_uppy_downy_matrix;
		}
		// if (prev == 542061 && len == (t - 1)) {
		//   std::cout << "Bhere I am, t = " << t << std::endl;
		//   std::cout << "uppy_downy_count = " << uppy_downy_count << std::endl;
		//   std::cout << "i0_cur = " << i0_cur << ", cur = " << cur << std::endl;
		//   std::cout << "i0_prev = " << i0_prev << ", prev = " << prev << std::endl;
		// }
		uppy_downy_matrix(uppy_downy_count, 0) = prev;
		uppy_downy_matrix(uppy_downy_count, 1) = t - 1;
		uppy_downy_matrix(uppy_downy_count, 2) = len;
		uppy_downy_count++;
	    }
	  }
	  i0_prev++;
	}
	if (verbose) {
	  std::cout << "fill in" << std::endl;	//## now fill in what was not set
	}
	while((i0_cur <= (mspbwtL - 1))) {
	  if ((fg + i0_cur + 1) <= (K - 1)) {
	    cur = a(fg + i0_cur + 1, t); // ## go up, so subtract
	    ud_down_cur(i0_cur) = cur;
	    e1 = t;
	    if (!use_cols_to_use0) {
	      if (!use_XR) {
		while((1 <= e1) && (X(cur, e1 - 1) == Z(e1 - 1))) {
		  e1--;
		}
	      } else {
		while((1 <= e1) && (XR(cur, e1 - 1) == Z(e1 - 1))) {
		  e1--;
		}
	      }
	    } else {
	      if (!use_XR) {
		while((1 <= e1) && (X(cur, cols_to_use0(e1 - 1)) == Z(e1 - 1))) {
		  e1--;
		}
	      } else {
		while((1 <= e1) && (XR(cur, cols_to_use0(e1 - 1)) == Z(e1 - 1))) {
		  e1--;
		}
	      }
	    }
	    //## go backwards, sort out start, how many before
	    ud_down_length_cur(i0_cur) = t - e1;
	  } else {
	    ud_down_cur(i0_cur) = -1;
	    ud_down_length_cur(i0_cur) = -1;
	  }
	  i0_cur++;
	}
	ud_down_prev = ud_down_cur;
	ud_down_length_prev = ud_down_length_cur;
	// if (t >= 20 && t <= 25) {
	//   std::cout << "first 5 going up at this location" << std::endl;
	//   std::cout << ud_up_cur(0) << ", " << ud_up_cur(1) << ", " << ud_up_cur(2) << ", " << ud_up_cur(3) << ", " << ud_up_cur(4) << std::endl;
	//   std::cout << "and their lengths" << std::endl;
	//   std::cout << ud_up_length_cur(0) << ", " << ud_up_length_cur(1) << ", " << ud_up_length_cur(2) << ", " << ud_up_length_cur(3) << ", " << ud_up_length_cur(4) << std::endl;
	//   std::cout << "first 5 going down at this location" << std::endl;
	//   std::cout << ud_down_cur(0) << ", " << ud_down_cur(1) << ", " << ud_down_cur(2) << ", " << ud_down_cur(3) << ", " << ud_down_cur(4) << std::endl;
	//   std::cout << "and their lengths" << std::endl;
	//   std::cout << ud_down_length_cur(0) << ", " << ud_down_length_cur(1) << ", " << ud_down_length_cur(2) << ", " << ud_down_length_cur(3) << ", " << ud_down_length_cur(4) << std::endl;
	// }
	if (verbose) {
	  std::cout << "end of up and down scna" << std::endl;
	}
      } // end of do up and down scan
      //
      //
      //
      if (!(g1 > f1)) {
	if (verbose) {
	  std::cout << "save and restart" << std::endl;
	}
	//
	// save and re-start! first, save
	//
	if (verbose) {
	  std::cout << "save first, fc = " << fc << ", gc = " << gc << std::endl;
	}
	if (!do_up_and_down_scan) {	
	  for(k = fc; k <= (gc - 1); k++) {
	    if ((t - ec) > min_length) {
	      Rcpp::IntegerVector temp_vector(4);
	      temp_vector(0) = k;
	      temp_vector(1) = a(k, t - 1);
	      temp_vector(2) = ec + 1;
	      temp_vector(3) = t - 1;
	      top_matches_list.push_back(temp_vector);
	    }
	  }
	}
	//
	if (verbose) {
	  std::cout << "restart" << std::endl;
	}
	// now, re-start
	//
	e1 = d(f1, t) - 1;
	fc = f1;
	gc = g1;
	matches_lower = false;
	matches_upper = false;
	// std::cout << "e1=" << e1 << ", t = " << t << ", f1 = " << f1 << ", K = " << K << std::endl;	    	    	  	  	
	if ((e1 == t) && (f1 == K)) {
	  e1 = t - 1;
	}
	//
	// see R code for explanation / comments
	if (verbose) {
	  std::cout << "both" << std::endl;
	  std::cout << "f1 = " << f1 << ", e1 = " << e1 << std::endl;
	}
	//
	while((!matches_lower) && (!matches_upper)) {
	  if (f1 > 0) {
	    if (!use_cols_to_use0) {
              e1_local = e1;
	    } else {
              e1_local=cols_to_use0(e1);
	    }
	    if (use_XR) {
	      matches_upper = Z(e1) == XR(a(f1 - 1, t), e1_local);
	    } else {
	      matches_upper = Z(e1) == X(a(f1 - 1, t), e1_local);
	    }
	  } else {
	    matches_upper = false;
	  }
	  if (f1 < K) {
	    if (!use_cols_to_use0) {
              e1_local = e1;
	    } else {
              e1_local=cols_to_use0(e1);
	    }
	    if (use_XR) {
	      matches_lower = Z(e1) == XR(a(f1, t), e1_local);
	    } else {
	      matches_lower = Z(e1) == X(a(f1, t), e1_local);
	    }
	  } else {
	    matches_lower = false;
	  }
	  if (!matches_lower & !matches_upper) {
	    e1++;
	  }
	}
	//
	//this CAN happen, if there is a symbol mis-match, and have to go forward
	if (verbose) {
	  std::cout << "matches upper" << std::endl;
	}
	//
	if (matches_upper) {
	    f1--;
	    index=a(f1, t);
            if (!use_cols_to_use0) {
	      if (use_XR) {
                while (Z(e1 - 1) == XR(index, e1 - 1)) {
		    e1--;
	        }
	      } else {
                while (Z(e1 - 1) == X(index, e1 - 1)) {
		    e1--;
	        }
	      }
	    } else {
	      if (use_XR) {
	        while (Z(e1 - 1) == XR(index, cols_to_use0(e1 - 1))) {
		    e1--;
	        }
	      } else {
	        while (Z(e1 - 1) == X(index, cols_to_use0(e1 - 1))) {
		    e1--;
	        }
	      }
	    }
	    while (d(f1, t) <= e1) {
	      f1--;
	    }
	}
	if (verbose) {
	  std::cout << "matches lower" << std::endl;
	}
	if (matches_lower) {
	    g1++;
	    index=a(f1, t);
	    if (!use_cols_to_use0) {
	      if (use_XR) {
		while (Z(e1 - 1) == XR(index, e1 - 1)) {
		  e1--;
		}
	      } else {
		while (Z(e1 - 1) == X(index, e1 - 1)) {
		  e1--;
		}
	      }
	    } else {
	      if (use_XR) {
		while (Z(e1 - 1) == XR(index, cols_to_use0(e1 - 1))) {
		  e1--;
		}
	      } else {
		while (Z(e1 - 1) == X(index, cols_to_use0(e1 - 1))) {
		  e1--;
		}
	      }
	    }
	    while ((g1 < K) && (d(g1, t) <= e1)) {
	      g1++;
	    }
	}
	ec = e1;
	if (verbose) {
	  std::cout << "end of save and restart" << std::endl;
	}
      }
      fc = f1;
      gc = g1;
      e(t - 1) = ec;
      f(t - 1) = fc;
      g(t - 1) = gc;
    }
    //
    // done normal now wrap up
    //
    if (verbose) {
      std::cout << "done normal now wrap up" << std::endl;
    }
    //t++; // no need to increment, for loop does that already
    if (!do_up_and_down_scan) {
      if (fc < gc) {
	if (verbose) {
	  std::cout << "fc = " << fc << ", gc = " << gc << std::endl;
	}
	for(k = fc; k <= (gc - 1); k++) {
	  if ((t - ec) > min_length) {	  
	    Rcpp::IntegerVector temp_vector(4);
	    temp_vector(0) = k;
	    temp_vector(1) = a(k, t - 1);
	    temp_vector(2) = ec + 1;
	    temp_vector(3) = t - 1;
	    top_matches_list.push_back(temp_vector);
	  }
	}
      }
      //
      // build final matrix
      //
      if (verbose) {
	std::cout << "build final matrix" << std::endl;
      }
      Rcpp::NumericMatrix top_matches(top_matches_list.length(), 4);
      for(k = 0; k < top_matches_list.length(); k++) {
	Rcpp::IntegerVector temp_vector = top_matches_list(k);
	for(int j = 0; j < 4; j++) {
	  top_matches(k, j) = temp_vector(j);
	}
      }
      colnames(top_matches) = Rcpp::CharacterVector({"k0", "indexB0", "start1", "end1"});
      return(top_matches);
    } else {
      //
      // report everything
      // up
      if (verbose) {
	std::cout << "start final summary of up and down scan" << std::endl;
      }
      i0_cur = 0;
      while((i0_cur <= (mspbwtL - 1))) {
	prev = ud_up_cur(i0_cur);
	len = ud_up_length_cur(i0_cur); 
	if (mspbwtM <= len) {
	  // Rcpp::IntegerVector temp_vector(3);			  
	  // temp_vector(0) = prev;
	  // temp_vector(1) = t - 1;
	  // temp_vector(2) = len;
	  if (uppy_downy_count >= uppy_downy_matrix.nrow()) {
	    // double the size!
	    Rcpp::IntegerMatrix new_uppy_downy_matrix(uppy_downy_matrix.nrow() * 2, 3);
	    for(int i = 0; i < uppy_downy_matrix.nrow(); i++) {
	      new_uppy_downy_matrix(i, _) = uppy_downy_matrix(i, _);
	    }
	    uppy_downy_matrix=new_uppy_downy_matrix;
	  }
	  uppy_downy_matrix(uppy_downy_count, 0) = prev;
	  uppy_downy_matrix(uppy_downy_count, 1) = t - 1;
	  uppy_downy_matrix(uppy_downy_count, 2) = len;
	  uppy_downy_count++;
	  //uppy_downy_reporter.push_back(temp_vector);
	}
	i0_cur++;
      }
      // down
      i0_cur = 0;
      while((i0_cur <= (mspbwtL - 1))) {
	// up
	prev = ud_down_cur(i0_cur);
	len = ud_down_length_cur(i0_cur);
	if (mspbwtM <= len) {
	  // Rcpp::IntegerVector temp_vector(3);			  
	  // temp_vector(0) = prev;
	  // temp_vector(1) = t - 1;
	  // temp_vector(2) = len;
	  if (uppy_downy_count >= uppy_downy_matrix.nrow()) {
	    // double the size!
	    Rcpp::IntegerMatrix new_uppy_downy_matrix(uppy_downy_matrix.nrow() * 2, 3);
	    for(int i = 0; i < uppy_downy_matrix.nrow(); i++) {
	      new_uppy_downy_matrix(i, _) = uppy_downy_matrix(i, _);
	    }
	    uppy_downy_matrix=new_uppy_downy_matrix;
	  }
	  uppy_downy_matrix(uppy_downy_count, 0) = prev;
	  uppy_downy_matrix(uppy_downy_count, 1) = t - 1;
	  uppy_downy_matrix(uppy_downy_count, 2) = len;
	  uppy_downy_count++;
	  //uppy_downy_reporter.push_back(temp_vector);
	}
	i0_cur++;
      }
      // now re-size again
      if (verbose) {
	std::cout << "resize" << std::endl;
      }
      Rcpp::NumericMatrix final_uppy_downy_matrix(uppy_downy_count, 4);
      for(k = 0; k < uppy_downy_count; k++) {
	final_uppy_downy_matrix(k, 0) = uppy_downy_matrix(k, 0);
	final_uppy_downy_matrix(k, 1) = uppy_downy_matrix(k, 1) - uppy_downy_matrix(k, 2) + 1;	
	final_uppy_downy_matrix(k, 2) = uppy_downy_matrix(k, 1);
	final_uppy_downy_matrix(k, 3) = uppy_downy_matrix(k, 2);	
      }
      colnames(final_uppy_downy_matrix) = Rcpp::CharacterVector({"index0", "start1", "end1", "len1"});
      if (verbose) {
	std::cout << "done" << std::endl;
      }
      return(final_uppy_downy_matrix);
    }
}









