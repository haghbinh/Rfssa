#include <RcppArmadillo.h>
// w-inner product
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
//'@importFrom Rcpp sourceCpp
// [[Rcpp::export]]
double winprod(arma::mat x, arma::mat y, arma::vec w,  arma::mat G)
{
  int N = x.n_cols;
  double s = 0;
  arma::vec xi(x.n_rows),yi(y.n_rows);
  for(int i=0; i<N; i++ )
  {
    xi = x.col(i);
    yi = y.col(i);
    s += arma::as_scalar(xi.t()*G*yi*w(i));
  }
  return(s);
}
