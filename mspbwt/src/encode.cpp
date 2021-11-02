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
    Rcpp::NumericMatrix out_mat(n_rows, col_names.length());
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
		(u(e0 + 1) - u(s0)) == 0 | \
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
