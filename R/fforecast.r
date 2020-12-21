#--------------------------------------------------------------
#' Functional Singular Spectrum Analysis Recurrent Forecasting and Vector Forecasting
#'
#' This function performs functional singular spectrum analysis (FSSA) recurrent forecasting (R-forecasting) or vector forecasting (V-forecasting) of functional time series (\code{\link{fts}}) as described in
#' Trinka et al. (2020)
#' @return A list of objects of class \code{\link{fts}} where each fts corresponds to a group
#' @param U an object of class \code{\link{fssa}} that holds the decomposition
#' @param groups a list of numeric vectors where each vector includes indices of elementary components of a group used for reconstruction and forecasting
#' @param h an integer that specifies the forecast horizon
#' @param method a character string specifying the type of forecasting to perform either 'recurrent' or 'vector'
#' @param tol a double specifying the amount of tolerated error in the approximation of the operator formed using a Neumann series leveraged in both forecasting algorithms
#' \dontrun{
#'# FSSA Forecasting
#'require(fda)
#'require(Rfssa)
#'data("Callcenter")
#'## Define functional objects
#'D <- matrix(sqrt(Callcenter$calls),nrow = 240)
#'N <- ncol(D)
#'time <- seq(ISOdate(1999,1,1), ISOdate(1999,12,31), by="day")
#'K <- nrow(D)
#'u <- seq(0,K,length.out=K)
#'d <- 23
#'basis <- create.bspline.basis(c(min(u),max(u)),d)
#'Ysmooth <- smooth.basis(u,D,basis)
#'## Define functional time series
#'Y <- fts(Ysmooth$fd,time = time)
#'plot(Y,xlab="Time (6 minutes aggregated)", ylab="Square root of call numbers",main = "Call Center Data")
#'## Perfrom FSSA decomposition
#'L <- 28
#'U <- fssa(Y,L)
#'groups <- list(1:7,1,2:3,4:5,6:7)
#'## Perfrom FSSA R-forecast and FSSA V-forecast
#'pr_R <- fforecast(U = U, groups = groups, h = 30, method = "recurrent", tol = 10^-3)
#'plot(pr_R[[1]],xlab="Time (6 minutes aggregated)", ylab="Square root of call numbers",main = "Recurrent Forecast Group 1")
#'plot(pr_R[[2]],xlab="Time (6 minutes aggregated)", ylab="Square root of call numbers",main = "Recurrent Forecast Group 2")
#'plot(pr_R[[3]],xlab="Time (6 minutes aggregated)", ylab="Square root of call numbers",main = "Recurrent Forecast Group 3")
#'plot(pr_R[[4]],xlab="Time (6 minutes aggregated)", ylab="Square root of call numbers",main = "Recurrent Forecast Group 4")
#'plot(pr_R[[5]],xlab="Time (6 minutes aggregated)", ylab="Square root of call numbers",main = "Recurrent Forecast Group 5")
#'
#'pr_V <- fforecast(U = U, groups = groups, h = 30, method = "vector", tol = 10^-3)
#'plot(pr_V[[1]],xlab="Time (6 minutes aggregated)", ylab="Square root of call numbers",main = "Vector Forecast Group 1")
#'plot(pr_V[[2]],xlab="Time (6 minutes aggregated)", ylab="Square root of call numbers",main = "Vector Forecast Group 2")
#'plot(pr_V[[3]],xlab="Time (6 minutes aggregated)", ylab="Square root of call numbers",main = "Vector Forecast Group 3")
#'plot(pr_V[[4]],xlab="Time (6 minutes aggregated)", ylab="Square root of call numbers",main = "Vector Forecast Group 4")
#'plot(pr_V[[5]],xlab="Time (6 minutes aggregated)", ylab="Square root of call numbers",main = "Vector Forecast Group 5")
#'
#' }
#' @export
fforecast <- function(U, groups=list(c(1)), h = 1, method = "recurrent",tol=10^-3){
  out=list()
  for(a in 1:length(groups)){
    g = groups[[a]]
    # Define prediction space
    basis <- U[[1]]$basis
    N <- U$N
    d <- nrow(U[[1]]$coefs)
    L <- ncol(U[[1]]$coefs)
    K <- N-L+1
    k <- length(g)
    D=matrix(data=NA,nrow=d,ncol=k)
    for(n in 1:k) D[,n]=U[[g[n]]]$coefs[,L];

    # Define Truncated Neumann Series
    NU_1 <- D%*%t(D)
    NU <- list(); NU[[1]]=diag(1,d);
    norm <- 1
    l <- 0
    while(norm>tol){
      l=l+1
      norm<-norm(NU[[l]])
      NU[[(l+1)]]=NU_1%*%NU[[l]]
    }
    Neu <- NU[[1]]
    for(n in 2:l) Neu=Neu+NU[[n]]


    if(method == "recurrent"){
      #FSSA R-forecasting

      # Reconstruct signal
      Q = Rfssa::freconstruct(U, groups = list(g))
      fssa_fore=matrix(data=0,nrow=d,ncol=h)
      for(m in 1:h){
        for(j in 1:(L-1)){
          E_j=matrix(data=NA,nrow=d,ncol=k)
          for(n in 1:k){
            E_j[,n]=U[[g[n]]]$coefs[,j]
          }
          A_j=Neu%*%D%*%t(E_j)
          fssa_fore[,m]=fssa_fore[,m]+A_j%*%as.matrix(Q[[1]][[1]]$coefs[,(N+j-L+m)])

        }
        Q[[1]][[1]]$coefs=cbind(Q[[1]][[1]]$coefs,fssa_fore[,m])
      }
      out[[a]]<-Rfssa::fts(fd(coef = fssa_fore,basisobj = basis))
    }else if(method == "vector"){
      # FSSA V-forecasting
      F_matrix=matrix(data=NA,nrow=((L-1)*d),ncol=k)
      for(n in 1:k){ F_matrix[,n]=matrix(data=U[[g[n]]]$coefs[,1:(L-1)],nrow = ((L-1)*d),ncol = 1)}
      P=F_matrix%*%t(F_matrix)+F_matrix%*%t(D)%*%Neu%*%D%*%t(F_matrix)
      S=array(data=0,dim=c(d,(K+h),L))
      for (j in 1L:length(g)){ S[,(1:K),] <- S[,(1:K),] + ufproj(U, g[j], d)}
      for(m in 1:h){
        obs=matrix(data=S[,(K+(m-1)),2:L],nrow = ((L-1)*d),ncol = 1)
        pr=P%*%obs
        pr=matrix(data=pr,nrow=d,ncol=(L-1))
        pr_1=matrix(data=0,nrow=d,ncol=1)
        for(j in 1:(L-1)){
          E_j=matrix(data=NA,nrow=d,ncol=k)
          for(n in 1:k){
            E_j[,n]=U[[g[n]]]$coefs[,j]
          }
          A_j=Neu%*%D%*%t(E_j)

          pr_1[,1]=pr_1[,1]+A_j%*%S[,(K+m-1),(j+1)]

        }
        pr=cbind(pr,pr_1)
        S[,(K+m),]=pr

      }
      S <- fH(S, d)
      predictions=S[,(K+1):(K+h),L]

      out[[a]]=Rfssa::fts(fd(coef = predictions, basisobj = basis))

    }
  }
  return(out)
}
