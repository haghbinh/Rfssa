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
#' ## Multivariate FSSA Example on Bivariate Satelite Image Data
#' library(fda)
#' library(Rfssa)
#' data(NDVI)
#' data(EVI)
#' d <- 11
#' basis <- create.bspline.basis(c(0,1),d)
#' u <- seq(0,1,length.out = 512)
#' y_1 <- smooth.basis(u,as.matrix(NDVI),basis)$fd
#' y_2 <- smooth.basis(u,as.matrix(EVI),basis)$fd
#' y=list(y_1,y_2)
#' Y=fts(y)
#' plot(Y)
#' fL=45
#' fU=fssa(Y,fL)
#' plot(fU,d=10,type='values')
#' plot(fU,d=10,type='paired')
#' plot(fU,d=10,type='lheats', var = 1)
#' plot(fU,d=10,type='lcurves',var = 1)
#' plot(fU,d=10,type='lheats', var = 2)
#' plot(fU,d=10,type='lcurves',var = 2)
#' plot(fU,d=10,type='wcor')
#' plot(fU,d=10,type='periodogram')
#' plot(fU,d=10,type='vectors')

#' frecon <- freconstruct(U = fU, group = list(c(1),c(2,3),c(4)))
#' plot(frecon[[1]],npts = 100,type = '3Dsurface',var=1)
#' plot(frecon[[2]],npts = 100,type = '3Dsurface',var=1)
#' plot(frecon[[3]],npts = 100,type = '3Dsurface',var=1)
#' plot(frecon[[1]],npts = 100,type = '3Dsurface',var=2)
#' plot(frecon[[2]],npts = 100,type = '3Dsurface',var=2)
#' plot(frecon[[3]],npts = 100,type = '3Dsurface',var=2)
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


