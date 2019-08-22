#' Addition of Functional Time Series
#'
#' + is an S3 method that lets you perform functional time series addition and scalar addition
#' @return An object of class fts.
#'
#' \item{out}{A functional time series}
#' @param Y1 A Univariate or Multivariate functional time series or scalar.
#' @param Y2 A Univariate or Multivariate functional time series or scalar.
#' @importFrom fda fd inprod eval.fd smooth.basis is.fd
#' @examples
#'
#' station_data = read.csv(file = "C:/Users/jorda/Dropbox/Jordan_Trinka_Research/MFSSA/Data/Pedestrian/Pedestrian_2015/January_2015.csv",header = TRUE)
#' library(fda)
#' library(Rfssa)
#' x=1:24
#' y=station_data[1:nrow(station_data),33:34]
#' nr=24
#' d=11
#' basis <- create.fourier.basis(c(0,1),d)
#' u=seq(0,1,length.out = 24)
#' Y=list()
#' for(j in 1:ncol(y)){
#'
#'  Y[[j]]=smooth.basis(argvals = u,y = sqrt(matrix(data=y[,j],nr=nr)),fdParobj = basis)$fd
#'
#'
#' }
#'
#' Y=fts(Y)
#'
#' Y1=fts(Y[[1]])
#'
#' Y2=fts(Y[[2]])
#'
#' plot(Y1 + Y2)
#' plot(2 + Y2)
#' plot(Y1 + 3)
#'
#'
#'
#' @useDynLib Rfssa
#' @export
'+.fts' <- function(Y1,Y2){
  if(class(Y1)=='fts' && class(Y2)=='fts'){
    p=Y1$p
    out = list()
    for(j in 1:p){

      basis_j=Y1[[j]]$basis

      out[[j]]=fd(Y1[[j]]$coefs+Y2[[j]]$coefs,basis_j)

      }

    out=fts(out)
    return(out)
  }else if(class(Y1)=='numeric' && class(Y2)=='fts'){
    p=Y2$p
    out = list()
    for(j in 1:p){

      basis_j=Y2[[j]]$basis

      out[[j]]=fd(Y1+Y2[[j]]$coefs,basis_j)

    }

    out=fts(out)
    return(out)


  }else if(class(Y1)=='fts' && class(Y2)=='numeric'){
    p=Y1$p
    out = list()
    for(j in 1:p){

      basis_j=Y1[[j]]$basis

      out[[j]]=fd(Y1[[j]]$coefs+Y2,basis_j)

    }

    out=fts(out)
    return(out)


  }

}

#' Subtraction of Functional Time Series
#'
#' - is an S3 method that lets you perform functional time series subtraction and scalar subtraction
#' @return An object of class fts.
#'
#' \item{out}{A functional time series}
#' @param Y1 A Univariate or Multivariate functional time series or scalar.
#' @param Y2 A Univariate or Multivariate functional time series or scalar.
#' @importFrom fda fd inprod eval.fd smooth.basis is.fd
#' @examples
#'
#' station_data = read.csv(file = "C:/Users/jorda/Dropbox/Jordan_Trinka_Research/MFSSA/Data/Pedestrian/Pedestrian_2015/January_2015.csv",header = TRUE)
#' library(fda)
#' library(Rfssa)
#' x=1:24
#' y=station_data[1:nrow(station_data),33:34]
#' nr=24
#' d=11
#' basis <- create.fourier.basis(c(0,1),d)
#' u=seq(0,1,length.out = 24)
#' Y=list()
#' for(j in 1:ncol(y)){
#'
#'  Y[[j]]=smooth.basis(argvals = u,y = sqrt(matrix(data=y[,j],nr=nr)),fdParobj = basis)$fd
#'
#'
#' }
#'
#' Y=fts(Y)
#'
#' Y1=fts(Y[[1]])
#'
#' Y2=fts(Y[[2]])
#'
#' plot(Y1 - Y2)
#'
#' plot(2 - Y2)
#'
#' plot(Y1 - 3)
#'
#'
#' @useDynLib Rfssa
#' @export
'-.fts' <- function(Y1,Y2){
  if(class(Y1) == 'fts' && class(Y2) == 'fts'){
    p=Y1$p
    out = list()
    for(j in 1:p){

      basis_j=Y1[[j]]$basis

      out[[j]]=fd(Y1[[j]]$coefs-Y2[[j]]$coefs,basis_j)

    }

    out=fts(out)
    return(out)
  }else if(class(Y1) == 'fts' && class(Y2) == 'numeric'){

    p=Y1$p
    out = list()
    for(j in 1:p){

      basis_j=Y1[[j]]$basis

      out[[j]]=fd(Y1[[j]]$coefs-Y2,basis_j)

    }

    out=fts(out)
    return(out)

  }else if(class(Y1) == 'numeric' && class(Y2) == 'fts'){

    p=Y2$p
    out = list()
    for(j in 1:p){

      basis_j=Y2[[j]]$basis

      out[[j]]=fd(Y1-Y2[[j]]$coefs,basis_j)

    }

    out=fts(out)
    return(out)

  }

}


#' Multiplication of Functional Time Series
#'
#' * is an S3 method that lets you multiply functional time series and perform scalar multiplication of functional time series
#' @return An object of class fts.
#'
#' \item{out}{A functional time series}
#' @param Y1 A Univariate or Multivariate functional time series or scalar
#' @param Y2 A Univariate or Multivariate functional time series or scalar.
#' @importFrom fda fd inprod eval.fd smooth.basis is.fd
#' @examples
#'
#' station_data = read.csv(file = "C:/Users/jorda/Dropbox/Jordan_Trinka_Research/MFSSA/Data/Pedestrian/Pedestrian_2015/January_2015.csv",header = TRUE)
#' library(fda)
#' library(Rfssa)
#' x=1:24
#' y=station_data[1:nrow(station_data),33:34]
#' nr=24
#' d=11
#' basis <- create.fourier.basis(c(0,1),d)
#' u=seq(0,1,length.out = 24)
#' Y=list()
#' for(j in 1:ncol(y)){
#'
#'  Y[[j]]=smooth.basis(argvals = u,y = sqrt(matrix(data=y[,j],nr=nr)),fdParobj = basis)$fd
#'
#'
#' }
#'
#' Y=fts(Y)
#'
#' Y1=fts(Y[[1]])
#'
#' Y2=fts(Y[[2]])
#'
#' plot(Y1 * Y2)
#'
#' plot(2 * Y2)
#'
#' plot(Y1 * 3)
#'
#'
#' @useDynLib Rfssa
#' @export
'*.fts' <- function(Y1,Y2){
  if(class(Y1) == 'fts' && class(Y2) == 'fts'){
    p=Y1$p
    out = list()
    for(j in 1:p){

      basis_j=Y1[[j]]$basis

      out[[j]]=fd(Y1[[j]]$coefs*Y2[[j]]$coefs,basis_j)

    }

    out=fts(out)
    return(out)
  }else if(class(Y1) == 'fts' && class(Y2) == 'numeric'){

    p=Y1$p
    out = list()
    for(j in 1:p){

      basis_j=Y1[[j]]$basis

      out[[j]]=fd(Y1[[j]]$coefs*Y2,basis_j)

    }

    out=fts(out)
    return(out)

  }else if(class(Y1) == 'numeric' && class(Y2) == 'fts'){

    p=Y2$p
    out = list()
    for(j in 1:p){

      basis_j=Y2[[j]]$basis

      out[[j]]=fd(Y1*Y2[[j]]$coefs,basis_j)

    }

    out=fts(out)
    return(out)

  }

}


#' Indexing into Functional Time Series
#'
#' [ is an S3 method that lets you index into a functional time series
#' @return An object of class fts.
#'
#' \item{out}{A functional time series}
#' @param Y A Univariate of Multivariate functional time series.
#' @importFrom fda fd inprod eval.fd smooth.basis is.fd
#' @examples
#'
#' station_data = read.csv(file = "C:/Users/jorda/Dropbox/Jordan_Trinka_Research/MFSSA/Data/Pedestrian/Pedestrian_2015/January_2015.csv",header = TRUE)
#' library(fda)
#' library(Rfssa)
#' x=1:24
#' y=station_data[1:nrow(station_data),33]
#' nr=24
#' d=11
#' basis <- create.fourier.basis(c(0,1),d)
#' u=seq(0,1,length.out = 24)
#' Y=list()
#' for(j in 1:ncol(y)){
#'
#'  Y[[j]]=smooth.basis(argvals = u,y = sqrt(matrix(data=y[,j],nr=nr)),fdParobj = basis)$fd
#'
#'
#' }
#'
#' Y=fts(Y)
#'
#'
#' plot(Y[3:5])
#'
#'
#'
#'
#'
#' @useDynLib Rfssa
#' @export
'[.fts'<-function(Y,i='index'){
  p=Y$p
  out=list()
  for(j in 1:p){

    basis_j=Y[[j]]$basis

    out[[j]]=fd(Y[[j]]$coefs[,i],basis_j)

  }

  out=fts(out)
  out$time=i
  return(out)

}


