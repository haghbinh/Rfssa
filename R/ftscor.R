#' Correlation for Functional Time Series Objects
#'
#' This function finds the correlation between univarite or multivariate functional time series objects
#' @return the correlation between functional time series.
#' @param Y1 a functional time series object of class \code{fts}.
#' @param Y2 a functional time series object of class \code{fts}.
#'
#'
#' @seealso \code{\link{fts}}
#' @export
cor.fts <- function(Y1, Y2) {
  N <- Y1$N
  w <- matrix(nrow=1,ncol=N,data=1)
  p = Y1$p
  G = list()
  Y1_list=list()
  Y2_list=list()
  for(j in 1:p){
    G[[j]] = inprod(Y1[[j]]$basis,Y1[[j]]$basis)
    Y1_list[[j]]=Y1[[j]]$coefs
    Y2_list[[j]]=Y2[[j]]$coefs
  }

  wcor <- mwinprod(Y1_list, Y2_list, w, G, p)/sqrt(mwinprod(Y1_list,Y1_list, w, G, p) * mwinprod(Y2_list,Y2_list, w, G, p))

  return(wcor)
}
