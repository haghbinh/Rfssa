#' Functional Singular Spectrum Analysis
#'
#' fssa is a function for decomposition stage (including embedding
#'  and  functional SVD steps) of a univariate or multivariate functional time series.
#' @return An object of class fssa, which is a list of
#' multivariate futional objects and the following components:
#' \item{values}{A numeric vector of eigenvalues.}
#' \item{L}{Window length.}
#' \item{N}{Length of time series.}
#' \item{Y}{The original functional time series.}
#' @param Y A functional time series.
#' @param L Window length.
#' @param type type of FSSA
#' @importFrom fda fd inprod eval.fd smooth.basis is.fd create.bspline.basis
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
#' Y <- fts(Ysmooth$fd)
#' ## fssa decomposition
#' L <- 28
#' U <- fssa(Y,L)
#' plot(U,d=13)
#' plot(U,d=9,type="lheats")
#' plot(U,d=9,type="lcurves")
#' plot(U,d=9,type="vectors")
#' plot(U,d=10,type="periodogram")
#' plot(U,d=10,type="paired")
#' plot(U,d=10,type="wcor")
#' ## fssa reconstruction
#' gr <- list(1,2:3,4:5,6:7,8:20)
#' Q <- freconstruct(U, gr)
#'
#'
#' layout(matrix(c(1,1,2,3,4,5,6,6),nr=2))
#' par(mar=c(2,1,2,2))
#' plot(Y,main="Call Numbers(Observed)")
#' plot(Q[[1]],main="1st Component")
#' plot(Q[[2]],main="2nd Component")
#' plot(Q[[3]],main="3rd Component")
#' plot(Q[[4]],main="4th Component")
#' plot(Q[[5]],main="5th Component(Noise)")
#'
#'
#' load('NDVINDVI.RData')
#' load('NDVIEVI.RData')

#'
#' d <- 11
#' basis <- create.bspline.basis(c(0,1),d)
#' u <- seq(0,1,length.out = 512)
#' y_1 <- smooth.basis(u,as.matrix(NDVI),basis)$fd
#' y_2 <- smooth.basis(u,as.matrix(EVI[,1:441]),basis)$fd
#' y=list(y_1,y_2)
#' library(Rfssa)
#' Y=fts(y)
#' plot(Y)
#' fL=45
#' fU=fssa(Y,fL)
#' plot(fU,d=10,type='values')
#' plot(fU,d=10,type='paired')
#' plot(fU,d=10,type='lheats', var = 1)
#' plot(fU,d=10,type='lheats', var = 2)
#' plot(fU,d=10,type='lcurves',var = 1)
#' plot(fU,d=10,type='lcurves',var = 2)
#' plot(fU,d=10,type='wcor')
#' plot(fU,d=10,type='periodogram')
#' plot(fU,d=10,type='vectors')

#' frecon <- freconstruct(U = fU, group = list(c(1),c(2,3),c(4)))
#' plot(frecon[[1]],npts = 100,type = '3Dsurface')
#' plot(frecon[[2]],npts = 100,type = '3Dsurface')
#' plot(frecon[[3]],npts = 100,type = '3Dsurface')
#'}
#' @useDynLib Rfssa
#' @export
fssa <- function(Y, L = NA) {
  if(is.fd(Y) & length(dim(Y$coefs)) == 2L )   Y <- fts(Y) else
    if(class(Y) != "fts") stop("The class of Y is not acceptable")
  if(is.na(L))  L <- floor(Y$N / 2L)
  if(Y$p==1){
    out <- ufssa(Y,L)
  } else if (Y$p > 1) {
    out <- mfssa(Y,L)
  } else stop("Dimension Error")
  class(out) <- "fssa"
  return(out)
}


