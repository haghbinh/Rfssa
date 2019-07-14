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
#' \item \code{"vectors"}  plot the eigenfunction's coefficients.(useful for the detection of period length).
#' \item \code{"meanvectors"}  plot the mean of eigenfunction's coefficients.(useful for the detection of period length).
#' \item \code{"meanpaired"} plot the pairs of mean of eigenfunction's coefficients. (useful for the detection of periodic components).
#' \item \code{"efunctions"} heatmap plot of eigenfunctions.(useful for the detection of period length).
#' \item \code{"efunctions2"} plot the eigenfunctions.(useful for the detection of meaningful patterns).
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
#' plot(U,d=4,type="efunctions")
#' plot(U,d=4,type="vectors")
#' plot(U,d=5,type="paired")
#' plot(U,d=5,type="wcor")
#' plot(U,d=5,type="meanvectors")
#' plot(U,d=5,type="efunctions2")
#' plot(U,d=5,type="meanpaired")
#' @seealso \code{\link{fssa}}, \code{\link{ftsplot}}
#' @export
plot.fssa <- function(x, d = length(x$values),
                      type = "values",...) {
  val <- sqrt(x$values)[1L:d]
  A <- val/sum(val)
  pr = round(A * 100L, 2L)
  main1 = paste0(1L:d, "(", pr,"%)")
  basis <- x[[1L]]$basis
  N <- x$N
  L <- x$L
  if (type == "values") {
    graphics::plot(val, type = "o", lwd = 2L,
         col = "dodgerblue3", pch = 19L,
         cex = 0.8, main = "Singular Values",
         ylab = " ", xlab = "Components")
  } else if (type == "paired") {
    n0 <- nrow(x$Y$coefs) * L
    x0 <- c(sapply(x[1L:d], function(x) as.vector(t(x$coefs))))
    D0 <- data.frame(x = x0[1L:((d - 1) * n0)], y = x0[(n0 + 1):(d * n0)])
    D0$group <- as.ordered(rep(paste(main1[1:(d - 1)], "vs", main1[2:d]), each = n0))
    p1 <- lattice::xyplot(x ~ y | group,
                          data = D0, xlab = "",
                          ylab = "", main = "Pairs of eigenvectors",
                          scales = list(x = list(at = NULL),
                                        y = list(at = NULL)),
                          as.table = TRUE, type = "l")
    graphics::plot(p1)
  } else if (type == "wcor") {
    W = fwcor(x, d)
    wplot(W)
  } else  if (type == "vectors") {
    n0 <- nrow(x$Y$coefs) * L
    x0 <- c(sapply(x[1L:d], function(x) as.vector(t(x$coefs))))
    D0 <- data.frame(x = x0,
                     time = rep(1L:n0, d))
    D0$group <- as.ordered(rep(main1,
                               each = n0))
    p1 <- lattice::xyplot(x ~ time |
                            group, data = D0, xlab = "",
                          ylab = "", main = "Eigenvectors",
                          scales = list(x = list(at = NULL),
                                        y = list(at = NULL)),
                          as.table = TRUE, type = "l")
    graphics::plot(p1)
  } else if (type == "meanvectors") {
    u <- basis$rangeval
    xindx <- seq(min(u), max(u),
                 length = 100)
    x0 <- c(sapply(x[1L:d],
                   function(x) colMeans(eval.fd(xindx,
                                                x))))
    D0 <- data.frame(x = x0,
                     time = rep(1L:L, d))
    D0$group <- as.ordered(rep(main1,
                               each = L))
    p1 <- lattice::xyplot(x ~ time |
                            group, data = D0, xlab = "",
                          ylab = "", main = "Meaned Eigenvectors",
                          scales = list(x = list(at = NULL),
                                        y = list(at = NULL)),
                          as.table = TRUE, type = "l")
    graphics::plot(p1)
  } else if (type == "meanpaired") {
    u <- basis$rangeval
    xindx <- seq(min(u), max(u),
                 length = 100L)
    x0 <- c(sapply(x[1L:d],
                   function(x) colMeans(eval.fd(xindx,
                                                x))))
    D0 <- data.frame(x = x0[1L:((d -
                                   1L) * L)], y = x0[(L +
                                                        1L):(d * L)])
    D0$group <- as.ordered(rep(paste(main1[1:(d -
                                                1L)], "vs", main1[2L:d]),
                               each = L))
    p1 <- lattice::xyplot(x ~ y | group,
                          data = D0, xlab = "",
                          ylab = "", main = "Meaned Paired Eigenvectors",
                          scales = list(x = list(at = NULL),
                                        y = list(at = NULL)),
                          as.table = TRUE, type = "l")
    graphics::plot(p1)
  }  else  if (type == "efunctions") {
    u <- basis$rangeval
    xindx <- seq(min(u), max(u),
                 length = 100L)
    n <- length(xindx)
    z0 <- lapply(x[1L:d], function(x) t(eval.fd(xindx,
                                                x)))
    z <- c(sapply(z0, function(x) as.vector(x)))
    D0 <- expand.grid(x = 1L:L,
                      y = 1L:n, group = 1L:d)
    D0$z <- z
    D0$group <- as.ordered(rep(main1,
                               each = L * n))
    p1 <- lattice::levelplot(z ~ x *
                               y | group, data = D0,
                             colorkey = TRUE, cuts = 50L,
                             xlab = "", ylab = "",
                             scales = list(x = list(at = NULL),
                                           y = list(at = NULL)),
                             aspect = "xy", as.table = TRUE,
                             main = "Pattern of Eigenfunctions",
                             col.regions = grDevices::heat.colors(100))
    graphics::plot(p1)
  } else if (type == "efunctions2") {
    col2 <- grDevices::rainbow(L)
    d1 <- floor(sqrt(d))
    d2 <- ifelse(d1^2 < d,
                 d1 + 1L, d1)
    graphics::par(mfrow = c(d1, d2),
        mar = c(2, 2, 3, 1),oma=c(2,2,7,1),cex.main=1.6)
    for (i in 1:d) graphics::plot(x[[i]],
                        lty = 1, xlab = "",
                        main = main1[i], ylab = "",
                        lwd = 2, col = col2)
    graphics::title("Pattern of Eigenfunctions",outer = TRUE)
    graphics::par(mfrow = c(1, 1))
  } else if (type == "rightvectors"){
    n0 <- nrow(x$Y$coefs) * K
    x0 <- c(V(U, d, K, L, Y))
    D0 <- data.frame(x = x0,
                     time = rep(1L:K, d))
    D0$group <- as.ordered(rep(main1,
                               each = K))
    p1 <- lattice::xyplot(x ~ time |
                            group, data = D0, xlab = "",
                          ylab = "", main = "Right Eigenvectors",
                          scales = list(x = list(at = NULL),
                                        y = list(at = NULL)),
                          as.table = TRUE, type = "l")
    graphics::plot(p1)
  } else if (type == "rightpairs"){
    x0 <- c(V(U, d, K, L, Y))
    D0 <- data.frame(x = x0[1L:((d -
                                   1L) * K)], y = x0[(K +
                                                        1L):(d * K)])
    D0$group <- as.ordered(rep(paste(main1[1:(d -
                                                1L)], "vs", main1[2L:d]),
                               each = K))
    p1 <- lattice::xyplot(x ~ y | group,
                          data = D0, xlab = "",
                          ylab = "", main = "Paired Eigenvectors (Rigth)",
                          scales = list(x = list(at = NULL),
                                        y = list(at = NULL)),
                          as.table = TRUE, type = "l")
    graphics::plot(p1)
  } else {
    stop("Unsupported type of fssa plot!")
  }
}
