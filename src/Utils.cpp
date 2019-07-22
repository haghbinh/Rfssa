#include <Rcpp.h>
#include <list>
using namespace Rcpp;
// Module  function
//' @importFrom Rcpp sourceCpp
//' @useDynLib Rfssa, .registration = TRUE
NumericVector mod(int i, int L) {
  NumericVector out(2);
  if( i%L == 0L) {
    out[1]=(int)i/L;
  } else {
    out[1]= (int)i/L+1L;
  }
  out[0]=i-(out[1]-1)*L;
  return out;
}


// Cofmath values calculation
// Transforming the vector c=(c_1,...,c_Ld) to a d*L matrix of coeficientsJT.
//'@importFrom Rcpp sourceCpp
// [[Rcpp::export]]
NumericMatrix Cofmat (int d,int L, NumericVector cx)
{ NumericMatrix S(d,L);
  for(int j=0 ; j < d; j++){
    for(int i=0; i < L; i++) {
      S(j,i) = cx(j* L + i);
    }
  }
  return S;
}

// Sij value
//' @importFrom Rcpp sourceCpp

double Csij(int i, int j, int K, int L, NumericMatrix B){
  NumericVector a1, a2 ;
  int i1,i2,j1,j2;
  double out(0);
  a1=mod(i,L);
  a2=mod(j,L);
  i1=a1[0];
  i2=a2[0];
  j1=a1[1];
  j2=a2[1];
  for(int ii = 0 ; ii< (K); ii++) {
    out += B(i1+ii-1,(j1-1))*B(i2+ii-1,(j2-1));
  }
  return out;
}


// S values calculation
//'@importFrom Rcpp sourceCpp
// [[Rcpp::export]]
NumericMatrix SS (int K,int L, NumericMatrix B, int d)
{ NumericMatrix S(L*d,L*d);
  for(int i=0 ; i < (L*d); i++){
    for(int j=0; j < (L*d); j++) {
      S(j,i)=Csij(i+1,j+1,K,L,B);
    }
  }
  return S;
}


// Gram matrix calculation
//
//'@importFrom Rcpp sourceCpp
// [[Rcpp::export]]
NumericMatrix Gram (int K,int L, NumericMatrix A, int d)
{
  NumericMatrix G(L*d,L*d);
  for(int i=0 ; i < (L*d); i++){
    for(int j=0; j < (L*d); j++) {
      NumericVector a1, a2 ;
      int ri, rj, qi, qj ;
      a1=mod(i+1,L);
      a2=mod(j+1,L);
      ri=a1[0];
      rj=a2[0];
      qi=a1[1];
      qj=a2[1];
      if (ri == rj) {
        G(i,j)=A(qi-1,qj-1);
      } else{
        G(i,j)=0;
      }
    }
  }
  return G;
}



// multi Sij value. This function stays basically the same except now we use two different inner product matrices for each variabljte.
//' @importFrom Rcpp sourceCpp
double Csmik(int k, int i, int d_j_k, int d_j_i, int K, int L, NumericMatrix B_j_k, NumericMatrix B_j_i){
  NumericVector k_vec, i_vec;
  int r_k,q_k,r_i,q_i;
  double out(0);
  k_vec=mod((k-d_j_k),L);
  i_vec=mod((i-d_j_i),L);
  r_k=k_vec[0];
  q_k=k_vec[1];
  r_i=i_vec[0];
  q_i=i_vec[1];
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

		B_j_k=findsubmat(k, p, B, shifter);

		d_j_k=findshift(k, p, shifter);

	  for(int i=0; i<(L*d_tilde); i++){


				B_j_i=findsubmat(i, p, B, shifter);

				d_j_i=findshift(i, p, shifter);


				S(k,i)=Csmik((k+1), (i+1), d_j_k, d_j_i, K, L, B_j_k, B_j_i);

		}

	}
	return S;
  }

// multi Gram matrix calculation. Finds the multivariate gram matrix using sub-matrices since the entry is zero if j_k is not jt_i or r_k is not r_i.
//
//'@importFrom Rcpp sourceCpp
// [[Rcpp::export]]
NumericMatrix Gramm (int K,int L, int p, int d_tilde, std::vector<NumericMatrix> A, NumericMatrix shifter, NumericMatrix d)
{
  NumericMatrix G(L*d_tilde,L*d_tilde);
  NumericMatrix G_j;
 for(int j = 0; j < p; j++){

	  G_j = Gram(K, L, A[j], (d[j+1]/L));

	  for(int k = (shifter(0,(j+1))-1); k < (shifter(1,(j+1))); k++){
		  for(int i = (shifter(0,(j+1))-1); i < (shifter(1,(j+1))); i++){

			  G(k, i) = G_j((k-((shifter(0,(j+1)))-1)), i-((shifter(0,(j+1)))-1));

		  }
	  }

  }


  return G;
}


