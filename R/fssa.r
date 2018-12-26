#' Functional SSA
#'
#' fssa is a function for decomposition stage (including embeding
#'  and  SVD step) of a functional time series.
#' @return An object of class fssa, which is a named list of
#' multivariate functional objects and the following components:
#' \item{values}{A numeric vector of eigenvalues.}
#' \item{L}{Window length.}
#' \item{N}{Length of time series.}
#' \item{Y}{The original functional time series.}
#' @param Y A functional time series.
#' @param L Window length.
#' @importFrom fda fd inprod eval.fd smooth.basis
#' @examples
#' ## Call Center Data
#' data("Callcenter")
#' library(fda)
#' D <- matrix(sqrt(Callcenter$calls),nrow = 240)
#' N <- ncol(D)
#' time <- 1:N
#' K <- nrow(D)
#' u <- seq(0,K,length.out =K)
#' d <- 22 #Optimal Number of basises
#' basis <- create.bspline.basis(c(min(u),max(u)),d)
#' Ysmooth <- smooth.basis(u,D,basis)
<<<<<<< HEAD
#' Y <- Ysmooth$fd
#' ## fssa decomposition
#' L <- 28
=======
#' ## fssa decomposition
#' L <- 30
>>>>>>> 6fb841ed4c67fc07967f1181d83281fc88fd5144
#' U <- fssa(Y,L)
#' plot(U,d=25)
#' plot(U,d=13,type="functions")
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
#'
#' layout(matrix(c(1,1,2,3,4,5,6,6),nr=2))
#' par(mar=c(2,1,2,2))
#' ftsplot(u,time,Y,space = 0.2,type=3,ylab = "",xlab = "Day",main = "Call Numbers(Observed)")
#' ftsplot(u,time,Q[[1]],space = 0.2,type=3,ylab = "",xlab = "Day",main = "1st Component")
#' ftsplot(u,time,Q[[2]],space = 0.2,type=3,ylab = "",xlab = "Day",main = "2nd Component")
#' ftsplot(u,time,Q[[3]],space = 0.2,type=3,ylab = "",xlab = "Day",main = "3rd Component")
#' ftsplot(u,time,Q[[4]],space = 0.2,type=3,ylab = "",xlab = "Day",main = "4th Component")
#' ftsplot(u,time,Q[[5]],space = 0.2,type=3,ylab = "",xlab = "Day",main = "5th Component(Noise)")
#'
#' @useDynLib Rfssa
#' @export
fssa <- function(Y, L = floor(dim(Y$coefs)[2L] / 2L)) {
    N <- dim(Y$coefs)[2]
    basis <- Y$basis
    d <- basis$nbasis
    K <- N - L + 1L
    B <- inprod(Y, basis)
    A <- inprod(basis, basis)
    S0 <- SS(K, L, B, d)
    H <- solve(Gram(K, L, A, d))
    Q <- eigen(H %*% S0)
    Q$vectors <- Re(Q$vectors)
    out <- list(NA)
    d1 <- sum(Re(Q$values) > 0.001)
    for (i in 1L:(d1)) out[[i]] <- fd(Cofmat(d, L, Q$vectors[, i]), basis)
    out$values <- Re(Q$values[1L:d1])
    out$L <- L
    out$N <- N
    out$Y <- Y
    class(out) <- "fssa"
    return(out)
}
