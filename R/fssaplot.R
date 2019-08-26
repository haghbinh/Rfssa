#--------------------------------------------------------------
#' Plotting fssa Objects
#'
#'  Plotting method for objects inheriting from class \code{\link{fssa}}.
#' @param x a functional singular value decomposition object, time series objects, usually inheriting from class "fssa".
#' @param d an integer which is the number of elementary components in the plot.
#' @param type what type of plot should be drawn. Possible types are:
#' \itemize{
#' \item \code{"values"} plot the square-root of singular values (default).
#' \item \code{"paired"} plot the pairs of eigenfunction's coefficients. (useful for the detection of periodic components).
#' \item \code{"wcor"} plot the W-correlation matrix for the reconstructed objects.
#' \item \code{"vectors"} plot the eigenfunction's coefficients.(useful for the detection of period length).
#' \item \code{"lcurves"} plot of the eigenfunctions.(useful for the detection of period length).
#' \item \code{"lheats"} heatmap plot the eigenfunctions.(useful for the detection of meaningful patterns).
#' \item \code{"periodogram"} periodogram plot.(useful for the detecting the frequencies of oscillations in functional data).
#' }
#' @param ... Arguments to be passed to methods, such as graphical parameters.
#' @examples
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
#' plot(U)
#' plot(U,d=4,type="lcurves")
#' plot(U,d=4,type="vectors")
#' plot(U,d=5,type="paired")
#' plot(U,d=5,type="wcor")
#' plot(U,d=5,type="lheats")
#' plot(U,d=5,type="periodogram")
#' @seealso \code{\link{fssa}}, \code{\link{ftsplot}}
#' @export
plot.fssa <- function(x, d = length(x$values),
                      type = "values",var=1L,ylab=NA) {
  val <- sqrt(x$values)[1L:d]
  p <- x$Y$p
  A <- val/sum(val)
  pr = round(A * 100L, 2L)
  main1 = paste0(1L:d, "(", pr,"%)")
  N <- x$N
  L <- x$L
  K <- N-L+1L
  if (type %in% c("lheats","lcurves")) {
    u <- x$Y$rangeval
    xindx <- seq(min(u), max(u),length = 100L)
    z0 <- list()
    for (i in 1:d){
      if(is.fd(x[[i]]))  x[[i]] <- list(x[[i]])
      z0[[i]] <- t(eval.fd(xindx,x[[i]][[var]]) )
    }
  }
  if (type == "values") {
    graphics::plot(val, type = "o", lwd = 2L,
         col = "dodgerblue3", pch = 19L,
         cex = 0.8, main = "Singular Values",
         ylab = " ", xlab = "Components")
  } else if (type == "wcor") {
    if(is.fd(x[[1]])) W <- ufwcor(x, d) else W <- mfwcor(x, d)
    wplot(W)
  }  else  if (type == "lheats") {
    n <- length(xindx)
    z <- c(sapply(z0, function(x) as.vector(x)))
    D0 <- expand.grid(x = 1L:L,
                      y = 1L:n, group = 1L:d)
    D0$z <- z
    D0$group <- as.ordered(rep(main1,
                               each = L * n))
    title0 <- "Singular functions"
    if(p>1) title0 <- paste(title0,"of the variable",
                            ifelse(is.na(ylab),var,ylab))
    p1 <- lattice::levelplot(z ~ x *
                               y | group, data = D0,
                             colorkey = TRUE, cuts = 50L,
                             xlab = "", ylab = "",
                             scales = list(x = list(at = NULL),
                                           y = list(at = NULL)),
                             aspect = "xy", as.table = TRUE,
                             main = title0,
                             col.regions = grDevices::heat.colors(100))
    graphics::plot(p1)
  } else if (type == "lcurves") {
    col2 <- grDevices::rainbow(L)
    d1 <- floor(sqrt(d))
    d2 <- ceiling(d/d1)
    graphics::par(mfrow = c(d1, d2),
        mar = c(2, 2, 3, 1),oma=c(2,2,7,1),cex.main=1.6)
    title0 <- "Singular functions"
    if(p>1) title0 <- paste(title0,"of the variable",
                            ifelse(is.na(ylab),var,ylab))

    for (i in 1:d){
      graphics::plot(x[[i]][[var]],
                        lty = 1, xlab = "",ylim=range(z0),
                        main = main1[i], ylab = "",
                        lwd = 2, col = col2)
    graphics::title(title0,outer = TRUE)
    }
    graphics::par(mfrow = c(1, 1))
  } else if (type == "vectors"){
    x0 <- c(x$RVectrs[,1L:d])
    D0 <- data.frame(x = x0,
                     time = rep(1L:K, d))
    D0$group <- as.ordered(rep(main1,
                               each = K))
    p1 <- lattice::xyplot(x ~ time |
                            group, data = D0, xlab = "",
                          ylab = "", main = "Singular vectors",
                          scales = list(x = list(at = NULL),
                                        y = list(at = NULL)),
                          as.table = TRUE, type = "l")
    graphics::plot(p1)
  } else if (type == "paired"){
    x0 <- c(x$RVectrs[,1L:d])
    D0 <- data.frame(x = x0[1L:((d -
                                   1L) * K)], y = x0[(K +
                                                        1L):(d * K)])
    D0$group <- as.ordered(rep(paste(main1[1:(d -
                                                1L)], "vs", main1[2L:d]),
                               each = K))
    p1 <- lattice::xyplot(x ~ y | group,
                          data = D0, xlab = "",
                          ylab = "", main = "Paired Singular vectors (Right)",
                          scales = list(x = list(at = NULL),
                                        y = list(at = NULL)),
                          as.table = TRUE, type = "l")
    graphics::plot(p1)
  } else if (type == "periodogram"){
    ff <- function(x) {
      I <- abs(fft(x)/sqrt(K))^2
      P = (4/K) * I
      return(P[1:(floor(K/2) + 1)])
    }
    x0 <- c(apply(x$RVectrs[,1L:d],2,ff))
    D0 <- data.frame(x = x0,
                     time = rep((0:floor(K/2))/K, d))
    D0$group <- as.ordered(rep(main1,
                               each = (floor(K/2) + 1)))
    p1 <- lattice::xyplot(x ~ time |
                            group, data = D0, xlab = "",
                          ylab = "", main = "Periodogram of Singular vectors",
                          scales = list(y = list(at = NULL)),
                          as.table = TRUE, type = "l")
    graphics::plot(p1)
    }else {
    stop("Unsupported type of fssa plot!")
  }
}
