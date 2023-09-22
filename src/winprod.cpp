// Include necessary headers.
#include <RcppArmadillo.h>

// Load the RcppArmadillo library.
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Define a function called 'winprod'.
//'@importFrom Rcpp sourceCpp
 // [[Rcpp::export]]
 double winprod(arma::mat x, arma::mat y, arma::vec w, arma::mat G) {
   // Get the number of columns in matrices x and y (assumed to be the same).
   int N = x.n_cols;

   // Initialize a variable 's' to store the result of the weighted inner product.
   double s = 0;

   // Define vectors 'xi' and 'yi' to store columns of matrices x and y.
   arma::vec xi(x.n_rows), yi(y.n_rows);

   // Loop through each column (i) from 0 to N-1.
   for (int i = 0; i < N; i++) {
     // Extract the i-th column from matrices x and y.
     xi = x.col(i);
     yi = y.col(i);

     // Calculate the weighted inner product for the i-th column and add it to 's'.
     // The weighted inner product involves matrix G, vectors xi, yi, and weight w(i).
     s += arma::as_scalar(xi.t() * G * yi * w(i));
   }

   // Return the computed result 's'.
   return s;
 }
