#include <RcppEigen.h>  // Include the RcppEigen library.
#include <list>         // Include the list library.

using namespace Rcpp;  // Use the Rcpp namespace.

// Declare the dependency on RcppEigen.
// [[Rcpp::depends(RcppEigen)]]

// Define the 'mod' function.
// This function calculates the module of 'i' with respect to 'L'.
//' @importFrom Rcpp sourceCpp
 //' @useDynLib Rfssa, .registration = TRUE
 NumericVector mod(int i, int L) {
   NumericVector out(2);  // Create a NumericVector 'out' with two elements.

   // Check if 'i' is divisible by 'L'.
   if (i % L == 0L) {
     out[1] = (int)i / L;  // Calculate the quotient.
   } else {
     out[1] = (int)i / L + 1L;  // Calculate the quotient with rounding up.
   }

   // Calculate the remainder.
   out[0] = i - (out[1] - 1) * L;

   return out;  // Return the 'out' vector containing the remainder and quotient.
 }



#include <Rcpp.h>
using namespace Rcpp;

// Define the exported Rcpp function.
// [[Rcpp::export]]
NumericMatrix Cofmat(int d, int L, NumericVector cx) {
  NumericMatrix S(d, L); // Create a d x L matrix to store coefficients.

  // Loop through each dimension (j) and each coefficient (i).
  for (int j = 0; j < d; j++) {
    for (int i = 0; i < L; i++) {
      // Calculate the index in the input vector cx and assign it to the corresponding element in S.
      S(j, i) = cx(j * L + i);
    }
  }

  // Return the resulting d x L coefficient matrix.
  return S;
}


// Define the exported Rcpp function.
// [[Rcpp::export]]
SEXP CalculateInverse(const Eigen::Map<Eigen::MatrixXd> A) {
  // Perform LU decomposition on the input matrix A.
  Eigen::PartialPivLU<Eigen::MatrixXd> lu(A);

  // Return the inverse of the LU-decomposed matrix as an SEXP object.
  return Rcpp::wrap(lu.inverse());
}



// Include the necessary headers for Rcpp and Eigen.
// [[Rcpp::depends(RcppEigen)]]

// Define the exported Rcpp function.
// [[Rcpp::export]]
SEXP AtimesB(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B) {
  // Multiply matrix A by matrix B and store the result in matrix C.
  Eigen::MatrixXd C = A * B;

  // Wrap matrix C into an R object and return it.
  return Rcpp::wrap(C);
}





// Include the necessary headers.
#include <Rcpp.h>
using namespace Rcpp;

// Define the C++ function.
//' @importFrom Rcpp sourceCpp
 double Csij(int i, int j, int K, int L, NumericMatrix B) {
   NumericVector a1, a2 ;
   int i1,i2,j1,j2;
   double out(0);

   // Calculate the modulo of i and j with respect to L.
   a1 = mod(i, L);
   a2 = mod(j, L);

   // Extract the individual components from the vectors a1 and a2.
   i1 = a1[0];
   i2 = a2[0];
   j1 = a1[1];
   j2 = a2[1];

   // Loop through ii from 0 to (K-1).
   for(int ii = 0; ii < K; ii++) {
     // Calculate the product of corresponding elements from matrix B
     // and add it to the result variable 'out'.
     out += B(i1 + ii - 1, j1 - 1) * B(i2 + ii - 1, j2 - 1);
   }

   // Return the final result.
   return out;
 }



// Load the Rcpp library.
//'@importFrom Rcpp sourceCpp
 // [[Rcpp::export]]
 Rcpp::NumericMatrix SS(int K, int L, Rcpp::NumericMatrix B, int d) {
   // Create a NumericMatrix S with dimensions (L*d) x (L*d).
   Rcpp::NumericMatrix S(L * d, L * d);

   // Loop through rows and columns of S.
   for (int i = 0; i < (L * d); i++) {
     for (int j = 0; j < (L * d); j++) {
       // Calculate the value for S(j, i) using the Csij function.
       S(j, i) = Csij(i + 1, j + 1, K, L, B);
     }
   }

   // Return the computed NumericMatrix S.
   return S;
 }




// Load the Rcpp library.
//'@importFrom Rcpp sourceCpp
 // [[Rcpp::export]]
 Rcpp::NumericMatrix Gram(int K, int L, Rcpp::NumericMatrix A, int d) {
   // Create a NumericMatrix G with dimensions (L*d) x (L*d).
   Rcpp::NumericMatrix G(L * d, L * d);

   // Loop through rows and columns of G.
   for (int i = 0; i < (L * d); i++) {
     for (int j = 0; j < (L * d); j++) {
       Rcpp::NumericVector a1, a2;
       int ri, rj, qi, qj;
       a1 = mod(i + 1, L);
       a2 = mod(j + 1, L);
       ri = a1[0];
       rj = a2[0];
       qi = a1[1];
       qj = a2[1];
       if (ri == rj) {
         // If ri is equal to rj, set G(i, j) to the corresponding value in A.
         G(i, j) = A(qi - 1, qj - 1);
       } else {
         // If ri is not equal to rj, set G(i, j) to 0.
         G(i, j) = 0;
       }
     }
   }

   // Return the computed NumericMatrix G.
   return G;
 }


// multi Sij value. This function stays basically the same except now we use two different inner product matrices for each variabljte.
//' @importFrom Rcpp sourceCpp
double Csmik(int k, int i, int d_j_k, int d_j_i, int K, int L, NumericMatrix B_j_k, NumericMatrix B_j_i){
  NumericVector k_vec, i_vec;
  int r_k,q_k,r_i,q_i;
  double out(0);
  k_vec = mod((k-d_j_k),L);
  i_vec = mod((i-d_j_i),L);
  r_k = k_vec[0];
  q_k = k_vec[1];
  r_i = i_vec[0];
  q_i = i_vec[1];
	for(int m = 0 ; m < (K); m++) {
		out += B_j_k(r_k+m-1,(q_k-1))*B_j_i(r_i+m-1,(q_i-1));
	}
  return out;
}



//  Basis Matrix Indexing calculation. This finds the proper basis matrix to use in S_0 calculation of S_k,i
//
//'@importFrom Rcpp sourceCpp
NumericMatrix findsubmat(int k, int p, std::vector<NumericMatrix> B, NumericMatrix shifter)
{
  int low_range_t;
  int high_range_t;
  NumericMatrix sub_mat;

  for(int t = 1; t <=p; t++){

	  low_range_t = shifter(0,t);
	  high_range_t = shifter(1,t);

	  if(low_range_t <= (k+1) && (k+1) <= high_range_t ){

		  sub_mat = B[(t-1)];

	  }
  }
  return sub_mat;
}

//  Basis shift Indexing calculation. Finds the unique shift, j_k, for a basis element k
//
//'@importFrom Rcpp sourceCpp
int findshift(int k, int p, NumericMatrix shifter)
{
  int low_range_t;
  int high_range_t;
  int shift = 0;
  for(int t = 1; t <=p; t++){
	  low_range_t = shifter(0,t);
	  high_range_t = shifter(1,t);
	  if(low_range_t <= (k+1) && (k+1) <= high_range_t ){
		  shift = shifter(1,(t-1));
	  }
  }
  return shift;
}


// multi test matrix calculation. Calculates elements of S_0 using shift paramejter between variables
//
//'@importFrom Rcpp sourceCpp
// [[Rcpp::export]]
NumericMatrix SSM (int K, int L, int d_tilde, int p, std::vector<NumericMatrix> B, NumericMatrix shifter)
{
  int d_j_k;
  int d_j_i;
  NumericMatrix B_j_k;
  NumericMatrix B_j_i;
  NumericMatrix S(L*d_tilde,L*d_tilde);
  for(int k=0; k<(L*d_tilde); k++){
		B_j_k = findsubmat(k, p, B, shifter);
		d_j_k = findshift(k, p, shifter);
	  for(int i = 0; i<(L*d_tilde); i++){
				B_j_i = findsubmat(i, p, B, shifter);
				d_j_i = findshift(i, p, shifter);
				S(k,i) = Csmik((k+1), (i+1), d_j_k, d_j_i, K, L, B_j_k, B_j_i);
		}
	}
	return S;
  }




// Load the Rcpp library.
//'@importFrom Rcpp sourceCpp
 // [[Rcpp::export]]
 Rcpp::NumericMatrix Gramm(int K, int L, int p, int d_tilde, std::vector<Rcpp::NumericMatrix> A, Rcpp::NumericMatrix shifter, Rcpp::NumericMatrix d) {
   // Create a NumericMatrix G with dimensions (L * d_tilde) x (L * d_tilde).
   Rcpp::NumericMatrix G(L * d_tilde, L * d_tilde);
   Rcpp::NumericMatrix G_j;

   // Loop through each value of j.
   for (int j = 0; j < p; j++) {
     // Calculate the Gram matrix for the current A[j].
     G_j = Gram(K, L, A[j], (d[j + 1] / L));

     // Loop through k and i within the specified ranges.
     for (int k = (shifter(0, (j + 1)) - 1); k < (shifter(1, (j + 1))); k++) {
       for (int i = (shifter(0, (j + 1)) - 1); i < (shifter(1, (j + 1))); i++) {
         // Copy values from G_j to G at the specified positions.
         G(k, i) = G_j((k - ((shifter(0, (j + 1))) - 1)), i - ((shifter(0, (j + 1))) - 1));
       }
     }
   }

   // Return the computed NumericMatrix G.
   return G;
 }
