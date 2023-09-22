#include <RcppArmadillo.h>

// Load the RcppArmadillo library.
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Define the exported Rcpp function.
// [[Rcpp::export]]
double HLinprod(arma::mat x, arma::mat y, arma::mat G) {
  int N = x.n_cols;  // Get the number of columns in matrices x and y.
  double s = 0;  // Initialize the result variable.

  arma::vec xi(x.n_rows);  // Create a column vector xi.
  arma::vec yi(y.n_rows);  // Create a column vector yi.

  // Loop through each column of x and y.
  for (int i = 0; i < N; i++) {
    xi = x.col(i);  // Get the i-th column of matrix x.
    yi = y.col(i);  // Get the i-th column of matrix y.

    // Calculate the inner product of xi and yi and add it to the result.
    s += arma::as_scalar(xi.t() * yi);
  }

  // Return the final result.
  return s;
}
