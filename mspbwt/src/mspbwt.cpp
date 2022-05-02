#include <Rcpp.h>
using namespace Rcpp;



// likely to change this function substantially!

//' @export
// [[Rcpp::export]]
Rcpp::List Rcpp_one_move_forward_buildindices(
    Rcpp::NumericMatrix X1C,
    Rcpp::IntegerMatrix a,
    Rcpp::IntegerMatrix usg,
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
  Rcpp::List to_return;
  return(to_return);
}
