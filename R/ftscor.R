#' Correlation for Functional Time Series Objects
#'
#' This function finds the correlation between univarite or multivariate functional time series (\code{\link{fts}}) objects.
#' @return a scalar that is the correlation between \code{\link{fts}} objects
#' @param Y1 an object of class \code{\link{fts}}
#' @param Y2 an object of class \code{\link{fts}}
#'
#'
#' @seealso \code{\link{fts}}
#'
#' @examples
#'
#' \dontrun{
#' require(fda)
#' require(Rfssa)
#' d <- 11
#' basis <- create.bspline.basis(c(0,1),d)
#' u <- seq(0,1,length.out = 512)
#' y_1 <- smooth.basis(u,as.matrix(NDVI),basis)$fd
#' y_2 <- smooth.basis(u,as.matrix(EVI),basis)$fd
#' Y1=fts(y_1)
#' Y2=fts(y_2)
#' cor.fts(Y1,Y2)
#' }
#'
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
