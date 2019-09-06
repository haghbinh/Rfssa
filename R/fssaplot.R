#--------------------------------------------------------------
#' Plotting fssa Objects
#'
#'  Plotting method for objects inheriting from class \code{\link{fssa}}.
#' @param x a functional singular value decomposition object, time series objects, usually inheriting from class "fssa".
#' @param d an integer which is the number of elementary components in the plot.
#' @param idx A Vector of indices of eigen elements to plot
#' @param idy A second vector of indices of eigen elements to plot(for type="paired")
#' @param groups A list or vector of indices determines grouping used for the decomposition(for type="wcor")
#' @param contrib logical. If 'TRUE' (the default), the contribution of the component to the total variance is displayed.
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
#' @param var an integer Specify the variable number.
#' @param ylab The character vector of name of variables.
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
#' @seealso \code{\link{fssa}}, \code{\link{plot.fts}}
#' @export
plot.fssa <- function(x, d = length(x$values),
                      idx = 1:d, idy = idx+1, contrib = TRUE,
                      groups = as.list(1:d),
                      type = "values",var=1L,ylab=NA, ...) {
  p <- x$Y$p
  A <- ((x$values)/sum(x$values))[1L:d]
  pr <- round(A * 100L, 2L)
  idx <- sort(idx)
  idy <- sort(idy)
  if(max(idx) > d | min(idx)< 1) stop("The idx must be subset of 1:d.")
  d_idx <- length(idx)
  if(contrib){
    main1 <- paste0(idx, "(", pr[idx],"%)")
    main2 <- paste0(idy, "(", pr[idy],"%)")
  } else{
    main1 <- paste(idx)
    main2 <- paste(idy)
  }
  N <- x$N
  L <- x$L
  K <- N-L+1L
  if (type %in% c("lheats","lcurves")) {
    u <- x$Y$rangeval
    xindx <- seq(min(u), max(u),length = 100L)
    z0 <- list()
    for (i in 1:d_idx){
      if(is.fd(x[[idx[i]]]))  x[[idx[i]]] <- list(x[[idx[i]]])
      z0[[i]] <- t(eval.fd(xindx,x[[idx[i]]][[var]]) )
    }
  }
  if (type == "values") {
    val <- sqrt(x$values)[idx]
    graphics::plot(idx, val, type = "o", lwd = 2L,
         col = "dodgerblue3", pch = 19L,
         cex = 0.8, main = "Singular Values",
         ylab = "norms", xlab = "Components")
  } else if (type == "wcor") {
    W <- fwcor(x, groups)
    wplot(W)
  }  else  if (type == "lheats") {
    n <- length(xindx)
    z <- c(sapply(z0, function(x) as.vector(x)))
    D0 <- expand.grid(x = 1L:L,
                      y = 1L:n, groups = idx)
    D0$z <- z
    D0$groups <- factor(rep(main1,
                               each = L * n), levels = main1)
    title0 <- "Singular functions"
    if(p>1) title0 <- paste(title0,"of the variable",
                            ifelse(is.na(ylab),var,ylab))
    p1 <- lattice::levelplot(z ~ x *
                               y | groups, data = D0,
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
    d1 <- floor(sqrt(d_idx))
    d2 <- ceiling(d_idx/d1)
    graphics::par(mfrow = c(d1, d2),
        mar = c(2, 2, 3, 1),oma=c(2,2,7,1),cex.main=1.6)
    title0 <- "Singular functions"
    if(p>1) title0 <- paste(title0,"of the variable",
                            ifelse(is.na(ylab),var,ylab))

    for (i in 1:d_idx){
      graphics::plot(x[[idx[i]]][[var]],
                        lty = 1, xlab = "",ylim=range(z0),
                        main = main1[i], ylab = "",
                        lwd = 2, col = col2)
    graphics::title(title0,outer = TRUE)
    }
    graphics::par(mfrow = c(1, 1))
  } else if (type == "vectors"){
    x0 <- c(x$RVectrs[,idx])
    D0 <- data.frame(x = x0,
                     time = rep(1L:K, d_idx))
    D0$groups <- factor(rep(main1,
                               each = K), levels = main1)
    p1 <- lattice::xyplot(x ~ time |
                            groups, data = D0, xlab = "",
                          ylab = "", main = "Singular vectors",
                          scales = list(x = list(at = NULL),
                                        y = list(at = NULL, relation="free")),
                          as.table = TRUE, type = "l")
    graphics::plot(p1)
  } else if (type == "paired"){
    d_idy <- length(idy)
    if(d_idx != d_idy) stop("The length of idx and idy must be same")
    x0 <- c(x$RVectrs[,idx])
    y0 <- c(x$RVectrs[,idy])
    D0 <- data.frame(x = x0, y = y0)
    main3 <- paste(main1, "vs", main2)
    D0$groups <- factor(rep(main3, each = K), levels = main3)
    p1 <- lattice::xyplot(x ~ y | groups,
                          data = D0, xlab = "",
                          ylab = "", main = "Paired Singular vectors (Right)",
                          scales = list(x = list(at = NULL, relation="free"),
                                        y = list(at = NULL, relation="free")),
                          as.table = TRUE, type = "l")
    graphics::plot(p1)
  } else if (type == "periodogram"){
    ff <- function(x) {
      I <- abs(fft(x)/sqrt(K))^2
      P = (4/K) * I
      return(P[1:(floor(K/2) + 1)])
    }
    x0 <- c(apply(x$RVectrs[,idx],2,ff))
    D0 <- data.frame(x = x0,
                     time = rep((0:floor(K/2))/K, d_idx))
    D0$groups <- factor(rep(main1,
                               each = (floor(K/2) + 1)), levels = main1)
    p1 <- lattice::xyplot(x ~ time |
                            groups, data = D0, xlab = "",
                          ylab = "", main = "Periodogram of Singular vectors",
                          scales = list(y = list(at = NULL, relation="free")),
                          as.table = TRUE, type = "l")
    graphics::plot(p1)
    }else {
    stop("Unsupported type of fssa plot!")
  }
}
