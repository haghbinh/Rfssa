#include <RcppArmadillo.h>// [[Rcpp::export]]


// Inner product between two elements of H^pL.
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
//'@importFrom Rcpp sourceCpp
// [[Rcpp::export]]
double HpLinprod(std::vector<arma::mat> X, std::vector<arma::mat> Y, std::vector<arma::mat> G,int p)
{
  double s = 0;
  arma::mat x;
  arma::mat y;
  //arma::mat g;
  for(int j=0; j<p; j++){

	x=X[j];
	y=Y[j];
	//g=G[j];
	int N = x.n_cols;
	arma::vec xi(x.n_rows),yi(y.n_rows);
	for(int i=0; i<N; i++)
	{
		xi = x.col(i);
		yi = y.col(i);
		//s += arma::as_scalar(xi.t()*g*yi);
		s += arma::as_scalar(xi.t()*yi);
	}
  }
  return(s);
}
