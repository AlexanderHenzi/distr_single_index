#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double get_crps(
    NumericMatrix cdf,
    NumericVector d_y,
    IntegerVector p_y,
    List p_x) {
  int n = p_x.length();
  int m = cdf.ncol() - 1;
  double out = 0.0;
  IntegerVector p_x_tmp;
  int p_y_tmp, n_p_x;
  for (int i = 0; i < n; i++) {
    p_x_tmp = p_x[i];
    n_p_x = p_x_tmp.length();
    for (int j = 0; j < n_p_x; j++) {
      p_y_tmp = p_y[p_x_tmp[j] - 1] - 1;
      for (int k = 0; k < p_y_tmp; k++) {
        out = out + pow(cdf(i, k), 2.0) * d_y[k];
      }
      for (int k = p_y_tmp; k < m; k++) {
        out = out + pow(cdf(i, k) - 1.0, 2.0) * d_y[k];
      }
    }
  }
  return out;
}
