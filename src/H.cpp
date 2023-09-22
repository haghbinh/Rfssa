#include <RcppArmadillo.h>

// Load the RcppArmadillo library.
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Define the exported Rcpp function.
// [[Rcpp::export]]
arma::mat H(arma::mat A) {
  int m = A.n_rows;  // Get the number of rows in matrix A.
  int n = A.n_cols;  // Get the number of columns in matrix A.
  bool cont = false;  // Initialize a flag for matrix transposition.

  // Check if the number of columns is greater than the number of rows.
  if (n > m) {
    // Swap m and n to make sure m <= n.
    m = m + n;
    n = m - n;
    m = m - n;

    // Transpose matrix A.
    A = A.t();
    cont = true;  // Set the flag for matrix transposition.
  }

  // Perform operations for the upper triangular part of the matrix A.
  for (int k = 1; k <= n; k++) {
    double S = 0;
    for (int i = 1; i <= k; i++)
      S += A(i - 1, k - i);
    double xm = S / k;
    for (int i = 1; i <= k; i++)
      A(i - 1, k - i) = xm;
  }

  // Check if m > n and perform operations for the lower triangular part.
  if (m > n) {
    for (int k = (n + 1); k <= m; k++) {
      double S = 0;
      for (int i = (k - n + 1); i <= k; i++)
        S += A(i - 1, k - i);
      double xm = S / n;
      for (int i = (k - n + 1); i <= k; i++)
        A(i - 1, k - i) = xm;
    }
  }

  // Perform operations for the remaining part of the matrix.
  for (int k = m + 1; k <= (n + m - 1); k++) {
    double S = 0;
    for (int i = (k - n + 1); i <= m; i++)
      S += A(i - 1, k - i);
    double xm = S / (m + n - k);
    for (int i = (k - n + 1); i <= m; i++)
      A(i - 1, k - i) = xm;
  }

  // If matrix A was transposed earlier, transpose it back.
  if (cont)
    A = A.t();

  // Return the modified matrix A.
  return A;
}
