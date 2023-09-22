#include <RcppArmadillo.h>

// Load the RcppArmadillo library.
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Define the exported Rcpp function.
// [[Rcpp::export]]
double HpLinprod(std::vector<arma::mat> X, std::vector<arma::mat> Y, std::vector<arma::mat> G, int p) {
  double s = 0;  // Initialize the result variable.
  arma::mat x;   // Declare a matrix x.
  arma::mat y;   // Declare a matrix y.

  // Loop through each element (j) in X, Y, and G vectors.
  for (int j = 0; j < p; j++) {
    x = X[j];  // Get the j-th element from vector X.
    y = Y[j];  // Get the j-th element from vector Y.
    int N = x.n_cols;  // Get the number of columns in matrices x and y.
    arma::vec xi(x.n_rows);  // Create a column vector xi.
    arma::vec yi(y.n_rows);  // Create a column vector yi.

    // Loop through each column (i) of matrices x and y.
    for (int i = 0; i < N; i++) {
      xi = x.col(i);  // Get the i-th column of matrix x.
      yi = y.col(i);  // Get the i-th column of matrix y.

      // Calculate the inner product of xi and yi and add it to the result.
      s += arma::as_scalar(xi.t() * yi);
    }
  }

  // Return the final result.
  return s;
}
