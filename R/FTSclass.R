#' Functional time series class
#'
#' The function fts is used to create functional time-series objects.

#' @param Y  an object of class fd or a list of objects of class fd.
#' @param time the vector of times at which a time series was sampled.

#' @export
fts <- function(Y,time = NA){
  if(is.list(Y) & !is.fd(Y)){
    p <- length(Y)
    NN <- dd <- rmin <- rmax <- NA
    for(i in 1:p) {
      NN[i] <- dim(Y[[i]]$coefs)[2]
      dd[i] <- dim(Y[[i]]$coefs)[1]
      rmin[i] <- Y[[i]]$basis$rangeval[1]
      rmax[i] <- Y[[i]]$basis$rangeval[2]
    }
    if(sum(NN-NN[1])!=0) stop("The time lengths are not equal.") else{
      N <- NN[1]
    }
    if(sum(rmin-rmin[1])!=0 | sum(rmax-rmax[1])!=0) stop("The domains of the functions are not same.") else{
      rangeval <- c(rmin[1],rmax[1])
    }
  } else if(is.fd(Y)){
    p <- 1
    N <- dim(Y$coefs)[2]
    dd <- dim(Y$coefs)[1]
    rangeval <- c(Y$basis$rangeval[1] , Y$basis$rangeval[2])
    Y <- list(Y)
  } else stop("Y must be an object of class fd or a list of class fd")
   if(is.na(time)) time <- 1:N else if(length(time)!= N) stop("The time length is not equal to data length.")
    Y$p <- p
    Y$N <- N
    Y$d <- dd
    Y$rangeval <- rangeval
    Y$time <- time
    class(Y) <- "fts"
    return(Y)
}



