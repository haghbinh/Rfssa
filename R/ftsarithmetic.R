#' Addition of Functional Time Series
#'
#' A method that lets you perform functional time series (\code{\link{fts}}) addition and scalar addition.
#' @return an object of class \code{\link{fts}}.
#'
#' @param Y1 an object of class \code{\link{fts}} or scalar
#' @param Y2 an object of class \code{\link{fts}} or scalar
#' @examples
#'
#' \dontrun{
#' require(fda)
#' require(Rfssa)
#' data(Callcenter) # Read data
#' u=seq(0,1,length.out=240) # Define domain of functional data
#' d=12 # number of basis elements
#' basis=create.bspline.basis(rangeval = c(0,1),nbasis = d) # create basis object
#' smooth.calls=smooth.basis(u, matrix(nrow=240,ncol=365,Callcenter$calls), basis)
#' Y=fts(smooth.calls$fd) # create functional time series
#' plot(Y)
#' Yplus=Y+Y # add the functional time series to itself
#' plot(Yplus)
#' Yplus2=Y+2 # add 2 to every term in the functional time series
#' plot(Yplus2)
#' }
#'
#' @seealso \code{\link{fts}}
#' @useDynLib Rfssa
#' @export
'+.fts' <- function(Y1,Y2){
  if(class(Y1)=='fts' && class(Y2)=='fts'){
    if(Y1$N!=Y2$N){

      stop("Functional time series are of different length")

    }
    if(Y1$p!=Y2$p){

      stop("Functional time series have different number of covariates")

    }
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
#' A method that lets you perform functional time series (\code{\link{fts}}) subtraction and scalar subtraction.
#' @return an object of class \code{\link{fts}}.
#'
#'
#' @param Y1 an object of class \code{\link{fts}} or scalar
#' @param Y2 an object of class \code{\link{fts}} or scalar
#' @examples
#'
#' \dontrun{
#' require(fda)
#' require(Rfssa)
#' data(Callcenter) # Read data
#' u=seq(0,1,length.out=240) # Define domain of functional data
#' d=12 # number of basis elements
#' basis=create.bspline.basis(rangeval = c(0,1),nbasis = d) # create basis object
#' smooth.calls=smooth.basis(u, matrix(nrow=240,ncol=365,Callcenter$calls), basis)
#' Y=fts(smooth.calls$fd) # create functional time series
#' plot(Y)
#' Yminus=Y[4:8]-Y[14:18] # subtract elements of the functional time series from each other
#' plot(Yminus)
#' Yminus2=Y-2 # add 2 to every term in the functional time series
#' plot(Yminus2)
#' }
#'
#' @seealso \code{\link{fts}}
#' @useDynLib Rfssa
#' @export
'-.fts' <- function(Y1,Y2){
  if(class(Y1) == 'fts' && class(Y2) == 'fts'){
    if(Y1$N!=Y2$N){

      stop("Functional time series are of different length")

    }
    if(Y1$p!=Y2$p){

      stop("Functional time series have different number of covariates")

    }
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
#' A method that lets you multiply functional time series (\code{\link{fts}}) and perform scalar multiplication of functional time series.
#' @return an object of class \code{\link{fts}}
#'
#'
#' @param Y1 an object of class \code{\link{fts}} or scalar
#' @param Y2 an object of class \code{\link{fts}} or scalar
#' @examples
#'
#' \dontrun{
#' require(fda)
#' require(Rfssa)
#' data(Callcenter) # Read data
#' u=seq(0,1,length.out=240) # Define domain of functional data
#' d=12 # number of basis elements
#' basis=create.bspline.basis(rangeval = c(0,1),nbasis = d) # create basis object
#' smooth.calls=smooth.basis(u, matrix(nrow=240,ncol=365,Callcenter$calls), basis)
#' Y=fts(smooth.calls$fd) # create functional time series
#' plot(Y)
#' Ytimes=Y*Y # elementwise multiplication of the functional time series with itself
#' plot(Ytimes)
#' Ytimes2=2*Y # multiply 2 with every term in the functional time series
#' plot(Ytimes2)
#' }
#'
#' @seealso \code{\link{fts}}
#' @useDynLib Rfssa
#' @export
'*.fts' <- function(Y1,Y2){
  if(class(Y1) == 'fts' && class(Y2) == 'fts'){
    if(Y1$N!=Y2$N){

      stop("Functional time series are of different length")

    }
    if(Y1$p!=Y2$p){

      stop("Functional time series have different number of covariates")

    }
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
#' A method that lets you index into a functional time series (\code{\link{fts}}).
#' @return an object of class \code{\link{fts}}
#'
#'
#' @param Y an object of class \code{\link{fts}}
#' @param i index
#' @examples
#'
#' \dontrun{
#' require(fda)
#' require(Rfssa)
#' data(Callcenter) # Read data
#' u=seq(0,1,length.out=240) # Define domain of functional data
#' d=12 # number of basis elements
#' basis=create.bspline.basis(rangeval = c(0,1),nbasis = d) # create basis object
#' smooth.calls=smooth.basis(u, matrix(nrow=240,ncol=365,Callcenter$calls), basis)
#' Y=fts(smooth.calls$fd) # create functional time series
#' Yind=Y[4:8] # take only the 4th through 8th functions
#' plot(Yind)
#' Yminus=Y[4:8]-Y[14:18] # subtract functions from each other
#' plot(Yminus)
#' }
#'
#'
#'
#' @seealso \code{\link{fts}}
#' @useDynLib Rfssa
#' @note can use ':' as an operator to specify a range of indices
#' @export
'[.fts'<-function(Y,i='index'){
  p=Y$p
  out=list()
  for(j in 1:p){

    basis_j=Y[[j]]$basis

    out[[j]]=fd(Y[[j]]$coefs[,i],basis_j)

  }

  out=fts(out)
  out$time=Y$time[i]
  return(out)

}


