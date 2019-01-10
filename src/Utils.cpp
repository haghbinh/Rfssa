#include <Rcpp.h>
using namespace Rcpp;
// Module  function
//' @importFrom Rcpp sourceCpp
//' @useDynLib Rfssa, .registration = TRUE
NumericVector mod(int K, int W) {
  NumericVector out(2);
  if( K%W == 0L) {
    out[1]=(int)K/W;
  } else {
    out[1]= (int)K/W+1L;
  }
  out[0]=K-(out[1]-1)*W;
  return out;
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


// Cofmath values calculation
// Transforming the vector c=(c_1,...,c_Ld) to a d*L matrix of coeficients.
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
