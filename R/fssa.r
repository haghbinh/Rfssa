#' Functional Singular Spectrum Analysis
#'
#' fssa is a function which performs the decomposition (including embedding
#'  and  functional SVD steps) stage for univariate functional singular spectrum analysis (ufssa)
#'  or multivariate functional singular spectrum analysis (mfssa)
#'  depending on whether the supplied input is a univariate or
#'  multivariate functional time series object.
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
#'
#' \dontrun{
#' ## Univariate FSSA Example on Callcenter data
#' data("Callcenter")
#' require(fda)
#' require(Rfssa)
#' ## Define functional objects
#' D <- matrix(sqrt(Callcenter$calls),nrow = 240)
#' N <- ncol(D)
#' time <- 1:N
#' K <- nrow(D)
#' u <- seq(0,K,length.out =K)
#' d <- 22 #Optimal Number of basises
#' basis <- create.bspline.basis(c(min(u),max(u)),d)
#' Ysmooth <- smooth.basis(u,D,basis)
#' ## Define functional time series
#' Y <- fts(Ysmooth$fd)
#' ## Univariate functional singular spectrum analysis
#' L <- 28
#' U <- fssa(Y,L)
#' plot(U,d=13)
#' plot(U,d=9,type="lheats")
#' plot(U,d=9,type="lcurves")
#' plot(U,d=9,type="vectors")
#' plot(U,d=10,type="periodogram")
#' plot(U,d=10,type="paired")
#' plot(U,d=10,type="wcor")
#' gr <- list(1,2:3,4:5,6:7,8:20)
#' Q <- freconstruct(U, gr)
#' plot(Y,main="Call Numbers(Observed)")
#' plot(Q[[1]],main="1st Component")
#' plot(Q[[2]],main="2nd Component")
#' plot(Q[[3]],main="3rd Component")
#' plot(Q[[4]],main="4th Component")
#' plot(Q[[5]],main="5th Component(Noise)")
#'
#'
#' ## Multivariate FSSA Example on Bivariate Satelite Image Data
#' require(fda)
#' require(Rfssa)
#' ## Raw image data
#' NDVI=Jambi$NDVI
#' EVI=Jambi$EVI
#' ## Kernel density estimation of pixel intensity
#' D0_NDVI <- matrix(NA,nrow = 512, ncol = 448)
#' D0_EVI <- matrix(NA,nrow =512, ncol = 448)
#' for(i in 1:448){
#'   D0_NDVI[,i] <- density(NDVI[,,i],from=0,to=1)$y
#'   D0_EVI[,i] <- density(EVI[,,i],from=0,to=1)$y
#' }
#' ## Define functional objects
#' d <- 11
#' basis <- create.bspline.basis(c(0,1),d)
#' u <- seq(0,1,length.out = 512)
#' y_NDVI <- smooth.basis(u,as.matrix(D0_NDVI),basis)$fd
#' y_EVI <- smooth.basis(u,as.matrix(D0_EVI),basis)$fd
#' y=list(y_NDVI,y_EVI)
#' ## Define functional time series
#' Y=fts(y)
#' plot(Y)
#' L=45
#' ## Multivariate functional singular spectrum analysis
#' U=fssa(Y,L)
#' plot(U,d=10,type='values')
#' plot(U,d=10,type='paired')
#' plot(U,d=10,type='lheats', var = 1)
#' plot(U,d=10,type='lcurves',var = 1)
#' plot(U,d=10,type='lheats', var = 2)
#' plot(U,d=10,type='lcurves',var = 2)
#' plot(U,d=10,type='wcor')
#' plot(U,d=10,type='periodogram')
#' plot(U,d=10,type='vectors')
#' recon <- freconstruct(U = U, group = list(c(1),c(2,3),c(4)))
#' plot(recon[[1]],npts = 100,type = '3Dsurface',var=1)
#' plot(recon[[2]],npts = 100,type = '3Dsurface',var=1)
#' plot(recon[[3]],npts = 100,type = '3Dsurface',var=1)
#' plot(recon[[1]],npts = 100,type = '3Dsurface',var=2)
#' plot(recon[[2]],npts = 100,type = '3Dsurface',var=2)
#' plot(recon[[3]],npts = 100,type = '3Dsurface',var=2)
#'
#'}
#' @useDynLib Rfssa
#' @export
fssa <- function(Y, L = NA, type="fssa") {
  if(is.fd(Y) & length(dim(Y$coefs)) == 2L )   Y <- fts(Y) else
    if(class(Y) != "fts") stop("The class of Y is not acceptable")
  if(is.na(L))  L <- floor(Y$N / 2L)
  if(Y$p==1 && type=="fssa"){
    out <- ufssa(Y,L)
  } else if (Y$p > 1 || type=="mfssa") {
    out <- mfssa(Y,L)
  } else stop("Error in Type or Dimension")
  class(out) <- "fssa"
  return(out)
}


