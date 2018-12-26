#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::export]]

arma::mat H (arma::mat A){
  int m=A.n_rows;
  int n=A.n_cols;
  bool cont= false;
  if(n>m)
  {
    m=m+n;
    n=m-n;
    m=m-n;
    A= A.t();
    cont=true;
  }
  for (int k =1; k <= n; k++ )
  {
    double S = 0;
    for (int i = 1; i <= k; i++ )
      S += A(i-1,k-i);
    double xm=S/k;
    for (int i = 1; i <= k; i++ )
      A(i-1,k-i)=xm;
  }
  if(m>n) {
    for (int k =(n+1); k <= m; k++ )
    {
      double S=0;
      for (int i= (k-n+1);i<= k; i++ )
        S += A(i-1,k-i);
      double xm=S/n;
      for (int i= (k-n+1);i<= k; i++ )
        A(i-1,k-i)=xm;
    }
  }
  for (int k =m+1; k <= (n+m-1); k++ )
  {
    double S=0;
    for (int i = (k-n+1); i <= m; i++ )
      S += A(i-1,k-i);
    double xm=S/(m+n-k);
    for (int i = (k-n+1); i <= m; i++ )
      A(i-1,k-i)=xm;
  }
  if(cont)  A=A.t();
  return A;
}


