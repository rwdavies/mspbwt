#include <Rcpp.h>
using namespace Rcpp;



//' @export
// [[Rcpp::export]]
Rcpp::List Rcpp_encode_maximal_column_of_u(
                                           Rcpp::IntegerVector u,
                                           int egs,
                                           bool efficient = true
                                           ) {
  int u_len = u.length();
  //int n_rows = std::ceil(float()u_len / float()egs);
  int n_rows = (u_len + egs - 1) / egs;
  Rcpp::CharacterVector col_names;
  Rcpp::IntegerVector out_vec(u_len);
  int col_start1 = 0;
  int col_start0 = 1;    
  int col_vec_pos;
  int col_value;
  if (efficient) {
    col_names = Rcpp::CharacterVector::create("value", "vec_pos");
	col_value = 0;
    col_vec_pos = 1;
  } else {
    col_names = Rcpp::CharacterVector::create("start1", "start0", "value", "vec_pos");
    col_value = 2;
    col_vec_pos = 3;
  }
  //Rcpp::NumericMatrix out_mat(n_rows, col_names.length());
  Rcpp::IntegerMatrix out_mat(n_rows, col_names.length());    
  colnames(out_mat) = col_names;
  //
  //
  int vec_pos = 0;
  int i, s0, e0, e0c, cs, j, d;
  bool do_encoding_this_SNP, rt, finalb;
  //
  for(i = 0; i < n_rows; i++) {
    s0 = egs * i;
    e0 = egs * (i + 1) - 1;
	e0c = e0;
    if ((e0 + 1) > u_len) {
	  e0 = u_len - 1;
	  e0c = e0 - 1;
	}
	if (!efficient) {
	  out_mat(i, col_start1) =  s0 + 1; // 1-based start
	  out_mat(i, col_start0) =  s0; // 0-based start	  
	}
    out_mat(i, col_value) = u(s0);
	do_encoding_this_SNP = true;
	if (i < (n_rows - 1)) {
      if (
          (u(e0 + 1) - u(s0)) == 0 || \
          (u(e0 + 1) - u(s0)) == (e0 + 1 - s0)
          ) {
        out_mat(i, col_vec_pos) = vec_pos - 1;
        do_encoding_this_SNP = false;
      }
	}
	//std::cout << "i = " << i << ", do_encoding_this_SNP = " << do_encoding_this_SNP << std::endl;
	if (do_encoding_this_SNP) {
      // now build runs from this
      cs = 0;
      rt = true; // true = +1, FALSE = 0+
      for(j = s0; j <= e0; j++) {
		if (j == (u_len - 1)) {
		  d = 2; // trigger storage on "final" one
		} else {
          d = u(j + 1) - u(j);
		}
		finalb = j == e0;
        if (rt) {
		  if ((d == 1) & !finalb) {
            cs++;
          } else {
            out_vec(vec_pos) = cs;
            //names(out_vec)[vec_pos + 1] <- i                    
            vec_pos++;
            cs = 1;
            rt = false;
          }
        } else {
		  if ((d == 0) & !finalb) {
            cs++;
          } else {
            out_vec(vec_pos) = cs;
			//names(out_vec)[vec_pos + 1] <- i
			vec_pos++;
			cs = 1;
			rt = true;
          }
        }
      }
      out_mat(i, col_vec_pos) = vec_pos - 1;
	}
  }
  out_vec =  out_vec[Rcpp::Range(0, vec_pos - 1)];
  //
  return(Rcpp::List::create(
                            Rcpp::Named("out_mat") = out_mat,
                            Rcpp::Named("out_vec") = out_vec
                            ));
}




//' @export
// [[Rcpp::export]]
int Rcpp_decode_maximal_value_of_u(
                                   Rcpp::IntegerMatrix & out_mat,
                                   Rcpp::IntegerVector & out_vec,
                                   int v,
                                   int egs,
                                   bool do_checks = false
                                   ) {
  // remainder is how many more to go
  // e.g. 14 means you have to go 14 into the creation of this thing
  int i_row = (v) / egs;  
  int remainder  = v - (i_row) * egs;
  int out_mat_nrows = out_mat.nrow();
  if (out_mat.ncol() != 2) {
    std::cout << "Error in assumptions! Fix" << std::endl;
    return(-1);
  }
  int val, next_val, vec_pos, u, next_u;
  bool is_plus;
  if (remainder == 0) {
    return(out_mat(i_row, 0));
  } else {
    val = out_mat(i_row, 0);
    if ((i_row + 1) < out_mat_nrows) {
      // these are constant - so don't need to check in between!
      next_val = out_mat(i_row + 1, 0); // "value" - assume efficient version!
      if (next_val == val) {
        return(val);
      } else if (((next_val) - val) == egs) {
        return(val + remainder);
      }
    }
    if (i_row == 0) {
	  vec_pos = 0;
    } else {
	  vec_pos = out_mat(i_row - 1, 1) + 1; // ## 0-based
    }
    u = 0;
    is_plus = true;
    while(u < remainder) {
      next_u = u + out_vec(vec_pos);
      if (remainder <= next_u) {
        if (is_plus) {
          val += remainder - u;
        }
        return(val);
      } else {
        u = next_u;
        if (is_plus) {
          val += out_vec(vec_pos);
          is_plus = false;
        } else {
          is_plus = true;
        }
        vec_pos++;
      }
    }
  }
  return(-1);
}




//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector Rcpp_encode_minimal_column_of_u(
                                                    Rcpp::IntegerVector & u
                                                    ) {
  int u_len = u.length();
  Rcpp::IntegerVector out(u(u_len - 1) - u(0));
  int j = 0;
  for(int i = 0; i < (u_len - 1); i++) {
    if ((u(i + 1) - u(i)) == 1) {
      out(j) = i + 1;
      j++;
    }
  }
  return(out);
}


//' @export
// [[Rcpp::export]]
int Rcpp_decode_minimal_value_of_u(
                                   Rcpp::IntegerVector & x,
                                   int v
                                   ) {
  // x is
  // v is 0-based on position
  int i = 0;
  int r = x(i);
  int x_len = x.length();
  while(v >= r) {
    i++;
    if ((i + 1) > x_len) {
      return(i);
	}
	r = x(i);
  }
  return(i);
}





// [[Rcpp::export]]
int is_list(SEXP x) {
  return Rf_isList(x) ? 0 : -1;
}



//' @export
// [[Rcpp::export]]
int Rcpp_decode_value_of_usge(
                              Rcpp::List usge,
                              int s,
                              int v,
                              int egs
                              ) {
  if (Rf_isNewList(usge(s - 1))) {
    Rcpp::List l = usge[s - 1];
    Rcpp::IntegerMatrix out_mat = l[0];
    Rcpp::IntegerVector out_vec = l[1];
    return( Rcpp_decode_maximal_value_of_u(out_mat, out_vec, v, egs));
  } else {
    Rcpp::IntegerVector x = usge[s - 1];
    return(Rcpp_decode_minimal_value_of_u(x, v));
  }
}


//' @export
// [[Rcpp::export]]
int Rcpp_decode_value_of_usge_v2(
                                 Rcpp::List& usge,
                                 int s,
                                 int v,
                                 int egs
                                 ) {
  if (Rf_isNewList(usge(s - 1))) {
    Rcpp::List l = usge[s - 1];
    Rcpp::IntegerMatrix out_mat = l[0];
    Rcpp::IntegerVector out_vec = l[1];
    return(Rcpp_decode_maximal_value_of_u(out_mat, out_vec, v, egs));
  } else {
    Rcpp::IntegerVector x = usge[s - 1];
    return(Rcpp_decode_minimal_value_of_u(x, v));
  }
}

