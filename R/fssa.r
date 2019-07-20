#' Functional Singular Spectrum Analysis
#'
#' fssa is a function for decomposition stage (including embedding
#'  and  functional SVD steps) of a functional time series.
#' @return An object of class fssa, which is a list of
#' multivariate functional objects and the following components:
#' \item{values}{A numeric vector of eigenvalues.}
#' \item{L}{Window length.}
#' \item{N}{Length of time series.}
#' \item{Y}{The original functional time series.}
#' @param Y A functional time series.
#' @param L Window length.
#' @importFrom fda fd inprod eval.fd smooth.basis
#' @examples
#' \dontshow{
#' require(Rfssa)
#' require(fda)
#' n <- 50 # Number of points in each function.
#' d <- 9
#' N <- 60
#' sigma <- 0.5
#' set.seed(110)
#' E <- matrix(rnorm(N*d,0,sigma/sqrt(d)),ncol = N, nrow = d)
#' basis <- create.fourier.basis(c(0, 1), d)
#' Eps <- fd(E,basis)
#' om1 <- 1/10
#' om2 <- 1/4
#' f0 <- function(tau, t) 2*exp(-tau*t/10)
#' f1 <- function(tau, t) 0.2*exp(-tau^3) * cos(2 * pi * t * om1)
#' f2 <- function(tau, t) -0.2*exp(-tau^2) * cos(2 * pi * t * om2)
#' tau <- seq(0, 1, length = n)
#' t <- 1:N
#' f0_mat <- outer(tau, t, FUN = f0)
#' f0_fd <- smooth.basis(tau, f0_mat, basis)$fd
#' f1_mat <- outer(tau, t, FUN = f1)
#' f1_fd <- smooth.basis(tau, f1_mat, basis)$fd
#' f2_mat <- outer(tau, t, FUN = f2)
#' f2_fd <- smooth.basis(tau, f2_mat, basis)$fd
#' Y_fd <- f0_fd+f1_fd+f2_fd
#' L <-10
#' U <- fssa(Y_fd,L)
#' gr <- as.list(1:4)
#' Q <- freconstruct(U, gr)
#' ftsplot(tau,t,Q[[1]],space = 0.2,type=3,xlab = "Day")
#' }
#' \dontrun{
#' ## Call Center Data
#' data("Callcenter")
#' require(fda)
#' require(Rfssa)
#' D <- matrix(sqrt(Callcenter$calls),nrow = 240)
#' N <- ncol(D)
#' time <- 1:N
#' K <- nrow(D)
#' u <- seq(0,K,length.out =K)
#' d <- 22 #Optimal Number of basises
#' basis <- create.bspline.basis(c(min(u),max(u)),d)
#' Ysmooth <- smooth.basis(u,D,basis)
#' Y <- Ysmooth$fd
#' ## fssa decomposition
#' L <- 28
#' U <- fssa(Y,L)
#' plot(U,d=13)
#' plot(U,d=9,type="efunctions")
#' plot(U,d=9,type="efunctions2")
#' plot(U,d=9,type="vectors")
#' plot(U,d=10,type="meanvectors")
#' plot(U,d=10,type="paired")
#' plot(U,d=10,type="meanpaired")
#' plot(U,d=10,type="wcor")
#' ## fssa reconstruction
#' gr <- list(1,2:3,4:5,6:7,8:20)
#' Q <- freconstruct(U, gr)
#'
#' cols3 <- rainbow(N)
#'
#' layout(matrix(c(1,1,2,3,4,5,6,6),nr=2))
#' par(mar=c(2,1,2,2))
#' plot(Y,lty=1,xlab="",main="Call Numbers(Observed)"
#'     ,ylab="",col=cols3)
#'     plot(Q[[1]],lty=1,xlab="",main="1st Component"
#'     ,ylab="",lwd=1,col=cols3)
#'     plot(Q[[2]],lty=1,xlab="",main="2nd Component"
#'     ,ylab="",lwd=1,col=cols3)
#'     plot(Q[[3]],lty=1,xlab="",main="3rd Component"
#'     ,ylab="",lwd=1,col=cols3)
#'     plot(Q[[4]],lty=1,xlab="",main="4th Component"
#'     ,ylab="",lwd=1,col=cols3)
#'     plot(Q[[5]],lty=1,xlab="",main="5th Component(Noise)"
#'     ,ylab="",lwd=1,col=cols3)
#'
#'
#' layout(matrix(c(1,1,2,3,4,5,6,6),nr=2))
#' par(mar=c(2,1,2,2))
#' ftsplot(u,time,Y,space = 0.2,type=3,ylab = "",xlab = "Day",main = "Call Numbers(Observed)")
#' ftsplot(u,time,Q[[1]],space = 0.2,type=3,ylab = "",xlab = "Day",main = "1st Component")
#' ftsplot(u,time,Q[[2]],space = 0.2,type=3,ylab = "",xlab = "Day",main = "2nd Component")
#' ftsplot(u,time,Q[[3]],space = 0.2,type=3,ylab = "",xlab = "Day",main = "3rd Component")
#' ftsplot(u,time,Q[[4]],space = 0.2,type=3,ylab = "",xlab = "Day",main = "4th Component")
#' ftsplot(u,time,Q[[5]],space = 0.2,type=3,ylab = "",xlab = "Day",main = "5th Component(Noise)")
#'}
#' @useDynLib Rfssa
#' @export
fssa <- function(Y, L = floor(dim(Y$coefs)[2L] / 2L, type="fssa")) {

  if(length(dim(Y$coefs))>2 & type == "fssa") stop("univariate fd FTS object
                                                    is needed for the type fssa.")

  if(type == "fssa"){
    out <- ufssa(Y,L)
  } else if (type == "mfssa") {
    out <- mfssa(Y,L)
  }
  out$type <- type
  class(out) <- "fssa"
  return(out)
}
