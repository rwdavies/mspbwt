#include <Rcpp.h>
using namespace Rcpp;



// likely to change this function substantially!

//' @export
// [[Rcpp::export]]
Rcpp::List Rcpp_ms_BuildIndices_Algorithm5(
    Rcpp::NumericMatrix X1C,
    Rcpp::List all_symbols,
    bool verbose = false,
    bool do_checks = false,
    int egs = 100,
    int n_min_symbols = 100
) {
    //
  //   Rcpp::IntegerVector n_symbols_per_grid(all_symbols.length());
  //   for(int i = 0; i < all_symbols.length(); i++) {
  //   //Rcpp::NumericMatrix x = all_symbols(i);
  //   //n_symbols_per_grid(i) = x.nrow();
  //   n_symbols_per_grid(i) = all_symbols[[i]](3, 1);
  // }
  //std::cout << n_symbols_per_grid << std::endl;
  Rcpp::List to_return;
  return(to_return);
}
