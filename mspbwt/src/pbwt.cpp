#include <Rcpp.h>
using namespace Rcpp;



//' @export
// [[Rcpp::export]]
Rcpp::List BuildIndices_Algorithm5_Rcpp(Rcpp::IntegerMatrix& X) {
  //
  const int K = X.nrow();
  const int T = X.ncol();
  int k, i, t, match_start, p, q;
  //
  IntegerMatrix a(K, T + 1);
  NumericMatrix d(K + 1, T + 1);
  NumericMatrix u(K + 1, T);
  NumericMatrix v(K + 1, T);
  for(k = 0; k < K; k++) {
    a(k, 0) = k;
    u(k, 0) = 0;
    v(k, 0) = 0;
    d(k, 0) = 0;
    d(k, 1) = 0;
  }
  NumericVector dtemp(K);
  NumericVector c(T + 1);
  NumericVector b(K);
  //
  // build first columns
  //
  int uu = 0;
  int vv = 0;
  for(k = 0; k < K; k++) {
    if (X(k, 0) == 0) {
      a(uu, 1) = k;
      uu++;
    } else {
      b(vv) = k;
      vv++;
    }
    u(k + 1, 0) = uu;
    v(k + 1, 0) = vv;
  }
  if (vv > 0) {
    for(i = 0; i <= (vv - 1); i++) {
      a(uu + i, 1) = b(i);
    }
  }
  d(uu, 1) = 1;
  c(0) = uu;
  c(T) = 0;
  //
  // add sentinels
  //
  for(t = 0; t <= T; t++) {
    d(0, t) = t + 1;
    d(K, t) = t + 1;
  }
  //
  //
  // t is now 0-based NOT 1-based
  // so - when it used as an index - don't need +1 change
  //
  for(t = 1; t < T; t++) {      
    uu = 0;
	vv = 0;
    p = t + 1; // ## 0-based SNP index + 1
    q = t + 1; // ## 0-based SNP index + 1
	for(k = 0; k < K; k++) {
	  match_start = d(k, t);
	  if (match_start > p) {
	    p = match_start;
	  }
	  if (match_start > q) {
	    q = match_start;
	  }
	  if (X(a(k, t), t) == 0) {
	    a(uu, t + 1) =  a(k, t);
	    d(uu, t + 1) = p;
	    uu++;
	    p = 0;
	  } else {
	    b(vv) = a(k, t);
	    dtemp(vv) = q;
	    vv++;
	    q = 0;
	  }
	  u(k + 1, t) = uu;
	  v(k + 1, t) = vv;
    }
    // 
    if (vv > 0) {
	  for(i = 0; i <= (vv - 1); i++) {
        a(uu + i, t + 1) = b(i);
        d(uu + i, t + 1) = dtemp(i);
	  }
	}
	c(t) = uu;
  }
  //
  return(Rcpp::List::create(
                            Rcpp::Named("a") = a,
                            Rcpp::Named("u") = u,
                            Rcpp::Named("v") = v,
                            Rcpp::Named("c") = c,
                            Rcpp::Named("d") = d
                            ));
}




//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix MatchZ_Algorithm5_Rcpp(
                                           Rcpp::IntegerMatrix& X,
                                           Rcpp::IntegerMatrix& a,
                                           Rcpp::NumericMatrix& u,
                                           Rcpp::NumericMatrix& v,
                                           Rcpp::NumericVector& c,
                                           Rcpp::NumericMatrix& d,
                                           Rcpp::NumericVector& Z,
                                           bool verbose = false,
                                           bool do_checks = false,
                                           int start_rows = 100
                                           ) {
  // misc init
  int k, index;
  int t = 0;
  int i_top_matches = 0;
  // start small, build out periodically if needed
  Rcpp::NumericMatrix top_matches(start_rows, 4);
  top_matches.fill(0);
  //
  const int K = X.nrow();
  const int T = X.ncol();
  //
  //
  Rcpp::NumericVector e(T);
  e.fill(0);
  Rcpp::NumericVector f(T);
  f.fill(0);
  Rcpp::NumericVector g(T);
  g.fill(0);
  if (Z(0) == 0) {
    f(0) = 0;
    g(0) = c(0);
    e(0) = 0;
  } else {
    f(0) = c(0);
    g(0) = K;
    e(0) = 0;
  }
  double fc = f(0);
  double gc = g(0);
  double ec = e(0);
  double e1 = -1;
  double f1 = -1;
  double g1 = -1;
  //
  // t is now 0-based NOT 1-based
  // so - when it used as an index - don't need +1 change
  //
  for(t = 1; t < T; t++) {
    if (Z(t) == 0) {
      f1 = u(fc, t);
	  g1 = u(gc, t);
    } else {
	  f1 = v(fc, t) + c(t);	    
	  g1 = v(gc, t) + c(t);
    }
    if (verbose) {
      std::cout << "Start of loop t = " << t;
      std::cout << ", fc = " << fc;
      std::cout << ", gc = " << gc;
      std::cout << ", ec = " << ec;
      std::cout << ", Z(t) = " << Z(t);
      std::cout << ", f1 = " << f1;
      std::cout << ", g1 = " << g1;
      std::cout << ", e1 = " << e1;			
      std::cout << std::endl;
      //message(paste0("Start of loop t=", t, ", fc = ", fc, ", gc = ", gc, ", ec = ", ec, ", Z[t] = ", Z[t],", f1=", f1, ", g1=", g1, ", e1 = ", e1))
    } 
    if (g1 > f1) {
      // nothing to do
    } else {
      if (verbose) {
        std::cout << "fc=" << fc << ", gc - 1 = " << gc - 1 << std::endl;
      }
      for(k = fc; k <= (gc - 1); k++) {
        if (verbose) {
          std::cout << "Save top match" << std::endl;
        }
        int current_rows = top_matches.nrow();
        if (i_top_matches >= current_rows) {
          Rcpp::NumericMatrix top_matches_temp(2 * current_rows, 4);
          for(int i_row = 0; i_row < current_rows; i_row++) {
            top_matches_temp.row(i_row) = top_matches.row(i_row);
          }
          top_matches = top_matches_temp;
        }
        top_matches(i_top_matches, 0) = k;
        top_matches(i_top_matches, 1) = a(k, t);
        top_matches(i_top_matches, 2) = ec + 1;
        top_matches(i_top_matches, 3) = t;
        i_top_matches++;
        if (verbose) {
          std::cout << "Done saving top match" << std::endl;
        }
      }
      //
      //
      //
      e1 = d(f1, t + 1) - 1;
      if (((Z(e1) == 0) && (f1 > 0)) || (f1 == K)) {
	    f1 = g1 - 1;
	    index = a(f1, t + 1);
	    while (Z(e1 - 1) == X(index, e1 - 1)) {
	      e1 = e1 - 1;
	    }
	    while (d(f1, t + 1) <= e1) {
	      f1 = f1 - 1;
	    }
      } else if (f1 < K) {
	    g1 = f1 + 1;
	    index = a(f1, t + 1);
	    while (Z(e1 - 1) == X(index, e1 - 1)) {
	      e1 = e1 - 1;
	    }
	    while ((g1 < K) && (d(g1, t + 1) <= e1)) {
	      g1 = g1 + 1;
	    }
      }
      ec = e1;
    }
    fc = f1;
    gc = g1;
    e(t) = ec;
    f(t) = fc;
    g(t) = gc;
  }
  //
  //
  //
  // t++;
  if (fc < gc) {
	for(k = fc; k <= (gc - 1); k++) {
	  int current_rows = top_matches.nrow();
	  if (i_top_matches >= current_rows) {
	    // re-allocate here
	    Rcpp::NumericMatrix top_matches_temp(2 * current_rows, 4);
	    for(int i_row = 0; i_row < current_rows; i_row++) {
	      top_matches_temp.row(i_row) = top_matches.row(i_row);
	    }
	    top_matches = top_matches_temp;
	  }
	  top_matches(i_top_matches, 0) = k;
	  top_matches(i_top_matches, 1) = a(k, t);
	  top_matches(i_top_matches, 2) = ec + 1;
	  top_matches(i_top_matches, 3) = t;
	  i_top_matches++;
	}
  }
  //
  // now shrink to size
  //
  Rcpp::NumericMatrix top_matches_temp(i_top_matches, 4);
  for(int i_row = 0; i_row < i_top_matches; i_row++) {
    top_matches_temp.row(i_row) = top_matches.row(i_row);
  }
  top_matches = top_matches_temp;
  colnames(top_matches) = CharacterVector::create("k0", "indexB0", "start1", "end1");
  //
  // done
  //
  return top_matches;
}
