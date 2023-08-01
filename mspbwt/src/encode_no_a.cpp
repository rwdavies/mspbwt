#include <Rcpp.h>
using namespace Rcpp;


int Rcpp_decode_value_of_usge(
    Rcpp::List usge,
    int s,
    int v,
    int egs
);



//' @export
// [[Rcpp::export]]
int Rcpp_get_k_given_encoded_u(
    int s,
    int v2,
    Rcpp::List& usge,
    Rcpp::IntegerMatrix& C,
    int K,
    int egs,
    bool verbose = false
) {
  //
  if (Rf_isNewList(usge(s - 1))) {  
      Rcpp::List l = usge[s - 1];
      Rcpp::IntegerMatrix out_mat = l[0];
      Rcpp::IntegerVector out_vec = l[1];
      if (verbose) {
	std::cout << "in slow one" << std::endl;
      }
      //
      // first, need the 0-based first position lower than it in out_mat
      // call it i_row
      //
      int k = -1;
      bool done = false;
      int out_mat_nrow = out_mat.nrow();
      //
      int c = 0;
      int up_i = 0;
      int up_val = 0;
      int down_i = out_mat_nrow - 1;
      int down_val = C(s - 1, 1);
      double frac = (double(v2) - double(up_val)) / (double(down_val) - double(up_val));
      if (frac > 1.001) {
	std::cout << "here comes an error!" << std::endl;
	std::cout << "s=" << s << ", v2 = " << v2 << std::endl;
	std::cout << "frac=" << frac <<  std::endl;	
	std::cout << "C" << std::endl;
	std::cout << C << std::endl;
	std::cout << "C(s - 1, 1) = " << C(s - 1, 1) << std::endl;            
	//std::cout << "C[s - 1, 1] = " << C[s - 1, 1] << std::endl;
	std::cout << "out_mat_nrow = " << out_mat_nrow << std::endl;
	Rcpp::stop("Something went wrong, with the choice of frac");
      }
      // for reasons I do not fully understand, round parentheses not workign properly
      //std::cout << "C(s, 1) = " << C(s, 1) << std::endl;            
      //std::cout << "C[s, 1] = " << C[s, 1] << std::endl;      
      //std::cout << "down_val = " << down_val << std::endl;
      //std::cout << "double(down_val) - double(up_val) = " << double(down_val) - double(up_val) << std::endl;
      int i_row = int(frac * double(K) /double(egs));
      i_row = std::min(out_mat_nrow - 2, i_row);      
      if (verbose) {
	std::cout << "C" << std::endl;
	std::cout << C << std::endl;
	std::cout << "s=" << s <<  std::endl;	
	std::cout << "C(s - 1, 1) = " << C(s - 1, 1) << std::endl;            
	//std::cout << "C[s - 1, 1] = " << C[s - 1, 1] << std::endl;      
	std::cout << "frac = " << frac << std::endl;
	std::cout << "at start, initial i_row=" << i_row << std::endl;
      }
      if (v2 >= out_mat(out_mat_nrow - 1, 0) ) {
	done = true;
	i_row = out_mat_nrow - 1;
      } 
      while(!done) {
	// not entirely sure what to do about last entry, we'll see
        // here let's make an estimate
        // check below then above
	// if neither we're good
	if ((i_row + 1) > out_mat_nrow) {
	  Rcpp::stop("Something went wrong, with the choice of i_row. Good luck");
	}
	if (v2 < out_mat(i_row, 0)) {
	  down_i = i_row;
	  down_val = out_mat(i_row, 0);
	  int i_row_proposed = int((double(v2) - double(up_val)) / (double(down_val) - double(up_val)) * (double(down_i) - double(up_i))) + up_i;
	  i_row = std::min(i_row - 1, i_row_proposed);
	} else if (out_mat(i_row + 1, 0) <= v2) {
	  up_i = i_row;
	  up_val = out_mat(i_row + 1, 0);
	  int i_row_proposed = int((double(v2) - double(up_val)) / (double(down_val) - double(up_val)) * (double(down_i) - double(up_i))) + double(up_i);
	  i_row = std::max(i_row + 1, i_row_proposed);
	} else {
	  done = true;
	}
	c++;
	if (c > 100000) {
	  Rcpp::stop("counter overflow problem in get_k_given_encoded_u");
	}
      }
      if (verbose) {
	std::cout << "i_row = " << i_row << std::endl;
      }
      //
      //
      //
      if (i_row < (out_mat_nrow - 1)) {
	if ((out_mat(i_row + 1, 0) - out_mat(i_row, 0)) == egs) {
	  k = egs * (i_row) + v2 - out_mat(i_row, 0);
	  return(k);
	}
      }
      //
      //
      //
      int vec_pos = -1;
      if (i_row == 0) {
	vec_pos = 0;
      } else {
	vec_pos = out_mat(i_row - 1, 1) + 1; // second column i.e. 1 here
      }
      int val = out_mat(i_row, 0);
      int remainder = v2 - val;
      int u = 0;
      int steps = 0;
      bool is_plus = true;
      done = false;
      while(!done) {
	if (verbose) {
	  std::cout << "before: vec_pos = " << vec_pos << ", u = " <<  u <<  ", steps = " << steps << std::endl;
	}
	if (is_plus) {
	  u += out_vec(vec_pos);
	}
	steps += out_vec(vec_pos);
	if (verbose) {
	  std::cout << "consider:vec_pos = " <<  vec_pos<< ", u = " <<  u <<  ", steps = " <<  steps << std::endl;
	}
	if (u > remainder) {
	  k = egs * i_row + steps - (u - remainder);
	  done = true;
	}
	vec_pos++;
	is_plus = !is_plus;
	if (!done && (vec_pos) > (out_mat(i_row, 1))) {
	  k = egs * (i_row + 1) - 1;
	  done = true;
	}
      }
      if (verbose) {
	std::cout << "steps = " << steps << std::endl;
      }
      return(k);
    } else {
    Rcpp::IntegerVector x = usge[s - 1];
      return(x(v2) - 1);
  }
}




//' @export
// [[Rcpp::export]]
int rcpp_go_backwards_one_step(
    int g,
    int v,
    Rcpp::IntegerMatrix& C,
    Rcpp::List& usge,
    int egs,
    int K,
    bool verbose = false
) {
  // get s (1-based)
  int s = 1;
  int count = C(0, 1);
  while(v > (count - 1)) {
    s++;
    count += C(s - 1, 1);
    if (verbose) {
      //std::cout << "s = " << s << ", C[s - 1, 1] = " << C[s - 1, 1] << ", C(s - 1, 1) = " << C(s - 1, 1) << ", v = " << v << ", count = " << count << std::endl;
      std::cout << "s = " << s << ", C(s - 1, 1) = " << C(s - 1, 1) << ", v = " << v << ", count = " << count << std::endl;
    }
  }
  // remember s is always 1-based
  int v2 = 0;
  if (s == 1) {
    v2 = v;
  } else {
    v2 = v - (count - C(s - 1, 1));
  }
  if (verbose) {
    std::cout << "s = " << s << ", v2 = " << v2 << std::endl;
  }
  int k = Rcpp_get_k_given_encoded_u(s, v2, usge, C, K, egs, verbose);
  if (verbose) {
    std::cout << "k = " << k << std::endl;
  }
  return(k);
}



//' @export
// [[Rcpp::export]]
int Rcpp_find_index_backward(
    int g_in,
    int v_in,
    Rcpp::List& all_symbols,
    Rcpp::List& usge_all,
    int egs,
    int K,
    Rcpp::List& list_of_columns_of_A,
    bool verbose = false,
    bool use_list_of_columns_of_A = false
) {
  int v = v_in;
  int k = -1;
  for(int g = g_in - 1; g >= (-1); g--) {
    //
    if (use_list_of_columns_of_A) {
      if (verbose) {
	std::cout << "check list_of_columns_of_A" << std::endl;
      }
      IntegerVector l = list_of_columns_of_A[g + 1 + 1];
      if (l.length() > 1) {
	if (verbose) {
	  std::cout << "return l[v]" << std::endl;
	}
	return(l[v]);
      }
    }
    //
    if (verbose) {
      std::cout << "g = " << g << ", v = " << v << std::endl;
    }
    Rcpp::List usge = usge_all[g + 1];
    Rcpp::IntegerMatrix C = all_symbols[g + 1];
    //std::cout << "C" << std::endl;
    //std::cout << C << std::endl;      
    k = rcpp_go_backwards_one_step(g, v, C, usge, egs, K, false);
    v = k;
  }
  // first col is 0:(K - 1) so this is OK!
  return(k);
}









//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector Rcpp_get_f_given_Z(
    Rcpp::IntegerVector & Z,    
    Rcpp::List & all_symbols,
    Rcpp::List & usge_all,
    int egs
) { 
  //
  //
  int s, k, u, c, s2;
  //
  //
  int nGrids = all_symbols.length();
  Rcpp::IntegerVector f(nGrids);
  //
  // init
  //
  Rcpp::IntegerMatrix C0 = all_symbols[0];
  int Z_1 = Z(0);
  if (Z_1 == 0) {
    Z_1 = C0.nrow();
  }
  if (Z_1 > 1) {
    for(s = 1; s <= (Z_1 - 1); s++) {
      f(0) += C0(s - 1, 1);
    }
  }
  //
  // update
  //
  for(int g = 1; g <= (nGrids - 1); g++) {
    Rcpp::IntegerMatrix C = all_symbols[g];
    Rcpp::List usge = usge_all[g];    
    s = Z(g);
    k = f(g - 1);
    u = Rcpp_decode_value_of_usge(usge, s, k, egs);
    c = 0;
    if (s > 1) {
      for(s2 = 2; s2 <= s; s2++) {
	c += C(s2 - 1 - 1, 1);
      }
    }
    u += c;
    f(g) = u;
  }
  return(f);
}



// note where f is     f <- get_f_given_Z(Z, all_symbols, usge_all, egs) in R


//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix Rcpp_find_good_matches_without_a(
    Rcpp::IntegerVector & Z,    
    Rcpp::List & all_symbols,
    Rcpp::List & usge_all,
    int egs,
    int pbwtL,
    int pbwtM,
    int K,
    Rcpp::RawMatrix & hapMatcherR,
    Rcpp::IntegerVector & which_snps_in_hapMatcherR,
    Rcpp::List &list_of_columns_of_A,
    bool use_list_of_columns_of_A,
    bool verbose = false
) {
  //
  //int a = -1;
  Rcpp::IntegerVector f = Rcpp_get_f_given_Z(Z, all_symbols, usge_all, egs);
  //
  int k, fc, g, c_up, c_down, l, v_up, cur_v, prev_k, v_down, c_mat, ic, len, g2, index, v, Zc;
  bool done;
  int n_c_mat = 0;
  //
  Rcpp::IntegerMatrix mat_up(pbwtL, 3); // v, k, l
  Rcpp::IntegerMatrix mat_down(pbwtL, 3);
  Rcpp::IntegerMatrix mat_up_prev(pbwtL, 3);
  Rcpp::IntegerMatrix mat_down_prev(pbwtL, 3);
  //  nam <- c("v", "s", "k", "l", "trueA")
  Rcpp::List to_out(f.length());
  Rcpp::LogicalVector to_out_used(f.length());
  for(g = f.length() - 1; g >= 0; g--) {
    if (verbose) {
      std::cout << "main loop, g = " << g << std::endl;
    }
    fc = f(g);
    Zc = Z(g);
    Rcpp::IntegerMatrix C = all_symbols[g];    
    Rcpp::List usge = usge_all[g];
    c_up = 0;
    c_down = 0;
    for(l = 0; l <= (pbwtL - 1); l++) {
        v_up = fc - l - 1;
	if (v_up >= 0 && v_up <= (K - 1)) {
	    k = rcpp_go_backwards_one_step(g + 1, v_up, C, usge, egs, K, false);
	    mat_up(l, 0) = v_up; // v
	    mat_up(l, 1) = k; // k
	    mat_up(l, 2) = 0; // l
	    cur_v = v_up;
	    prev_k = mat_up_prev(c_up, 1);
	    if (cur_v == prev_k) {
	      mat_up(l, 2) = mat_up_prev(c_up, 2) + 1; // l
	      c_up++;
	    }
	} else {
	  mat_up(l, 0) = -1;
	  mat_up(l, 1) = -1;
	  mat_up(l, 2) = -1;
	}
	v_down = fc + l;
	if (v_down >= 0 && v_down <= (K - 1)) {
	    k = rcpp_go_backwards_one_step(g + 1, v_down, C, usge, egs, K, false);
	    mat_down(l, 0) = v_down;
	    mat_down(l, 1) = k;
	    mat_down(l, 2) = 0;
	    cur_v = v_down;
	    prev_k = mat_down_prev(c_down, 1);
	    if (cur_v == prev_k) {
	      mat_down(l, 2) = mat_down_prev(c_down, 2) + 1;
	      c_down++;
	    }
	} else {
	  mat_down(l, 0) = -1;
	  mat_down(l, 1) = -1;
	  mat_down(l, 2) = -1;
	}
    }
    //
    // check re save condition
    //
    if ((g < (f.length() - 1) && (c_up < pbwtL || c_down < pbwtL)) || g == 0) {
      if (g == 0) {
	c_up = 0;
	c_down = 0;
      }
      Rcpp::IntegerMatrix mat_out(2 * pbwtL - c_up - c_down, 3);
      c_mat = 0;
      // up
      if (c_up < pbwtL) {
	for(ic = c_up; ic <= pbwtL - 1; ic++) {
	  if (mat_up_prev(ic, 0) >= 0) {
	    v = mat_up_prev(ic, 0);  // v ie 0-based col 0
	    k = mat_up_prev(ic, 1); // k ie 0-based col 1
	    len = mat_up_prev(ic, 2); // l i.e. 0-based col 2
	    // if (verbose) {
	    //   std::cout << "Up, Find index backward: g = " << g << " v_in = " << k << std::endl;
	    // }
	    index = Rcpp_find_index_backward(g, k, all_symbols, usge_all, egs, K, list_of_columns_of_A, false, use_list_of_columns_of_A);
	    g2 = g;
	    done = false;
	    while(!done & (g2 >= 0)) {
	      if (hapMatcherR(index, which_snps_in_hapMatcherR(g2) - 1) == Z[g2]) {
		len++;
		g2--;
	      } else {
		done = true;
	      }
	    }
	    if (len > pbwtM) {
	      mat_out(c_mat, 0) = g2 + 1;
	      mat_out(c_mat, 1) = index;
	      mat_out(c_mat, 2) = len;
	      c_mat++;
	    }
	  }
	}
      }
      // down
      if (c_down < pbwtL) {
	for(ic = c_down; ic <= pbwtL - 1; ic++) {
	  if (mat_down_prev(ic, 0) >= 0) {
	    v = mat_down_prev(ic, 0);  // v ie 0-based col 0
	    k = mat_down_prev(ic, 1); // k ie 0-based col 1
	    len = mat_down_prev(ic, 2); // l i.e. 0-based col 2
	    // if (verbose) {
	    //   std::cout << "Down, Find index backward: g = " << g << " v_in = " << k << std::endl;
	    // }
	    index = Rcpp_find_index_backward(g, k, all_symbols, usge_all, egs, K, list_of_columns_of_A, false, use_list_of_columns_of_A);
	    g2 = g;
	    done = false;
	    while(!done & (g2 >= 0)) {
	      if (hapMatcherR(index, which_snps_in_hapMatcherR(g2) - 1) == Z[g2]) {
		len++;
		g2--;
	      } else {
		done = true;
	      }
	    }
	    if (len > pbwtM) {
	      mat_out(c_mat, 0) = g2 + 1;
	      mat_out(c_mat, 1) = index;
	      mat_out(c_mat, 2) = len;
	      c_mat++;
	    }
	  }
	}
      }
      if (c_mat > 0) {
	n_c_mat = n_c_mat + c_mat;	
	Rcpp::IntegerMatrix mat_out2(c_mat, 3);
	for(ic = 0; ic < c_mat; ic++) {
	  mat_out2(ic, 0) = mat_out(ic, 0);
	  mat_out2(ic, 1) = mat_out(ic, 1);
	  mat_out2(ic, 2) = mat_out(ic, 2);	  
	}
	to_out[g] = mat_out2;
	to_out_used[g] = true;
      }
    }
    for(ic = 0; ic < pbwtL; ic++) {
      mat_up_prev(ic, 0) = mat_up(ic, 0);
      mat_up_prev(ic, 1) = mat_up(ic, 1);
      mat_up_prev(ic, 2) = mat_up(ic, 2);
      mat_down_prev(ic, 0) = mat_down(ic, 0);
      mat_down_prev(ic, 1) = mat_down(ic, 1);
      mat_down_prev(ic, 2) = mat_down(ic, 2);      
    }
  }
  //
  // merge at the end
  //
  if (verbose) {
    std::cout << "Merge at the end" << std::endl;
  }
  Rcpp::IntegerMatrix final_mat_out(n_c_mat, 3);
  colnames(final_mat_out) = Rcpp::CharacterVector({"g", "index", "len"});
  if (n_c_mat == 0) {
    return(final_mat_out);
  }
  c_mat = 0;
  for(g = 0; g < f.length(); g++) {
    if (to_out_used[g]) {
      Rcpp::IntegerMatrix m = to_out[g];
      for(ic = 0; ic < m.nrow(); ic++) {
	final_mat_out(c_mat + ic, 0) = m(ic, 0);
	final_mat_out(c_mat + ic, 1) = m(ic, 1);
	final_mat_out(c_mat + ic, 2) = m(ic, 2);	
      }
      c_mat += m.nrow();
    }
  }
  return(final_mat_out);
}
