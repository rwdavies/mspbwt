#include <Rcpp.h>
using namespace Rcpp;


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
	std::cout << "C[s - 1, 1] = " << C[s - 1, 1] << std::endl;
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
	std::cout << "C[s - 1, 1] = " << C[s - 1, 1] << std::endl;      
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
      std::cout << "s = " << s << ", C[s - 1, 1] = " << C[s - 1, 1] << ", C(s - 1, 1) = " << C(s - 1, 1) << ", v = " << v << ", count = " << count << std::endl;
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
      Rcpp::List l = list_of_columns_of_A[g + 1 + 1];
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
