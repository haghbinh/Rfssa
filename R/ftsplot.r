#--------------------------------------------------------------
#' Functional Time Series Plots
#'
#' Novel different type of plots to visualize the functional time series data objects.
#'
#' @param x a numeric vector of cordinate x.
#' @param y a numeric vector of cordinate y.
#' @param X a functional time series object of class \code{\link[fda]{fd}}.
#' @param type what type of plot should be drawn. Possible types are
#'     "type=1" for ribbon3D plot with curtain, "type=2" for ribbon3D plot
#'      without curtain, "type=3" for image2D plot.
#' @param zlab The lable of z axis.
#' @param xlab The lable of x axis.
#' @param ylab The lable of y axis.
#' @param space The amount of space (as a fraction of the average ribbon width) left between ribbons.
#' @param main The main title.
#' @examples
#' data("Callcenter")
#' library(fda)
#' D <- matrix(sqrt(Callcenter$calls),nrow = 240)
#' N <- ncol(D)
#' time <- 1:30
#' K <- nrow(D)
#' u <- seq(0,K,length.out =K)
#' d <- 22 #Optimal Number of basises
#' basis <- create.bspline.basis(c(min(u),max(u)),d)
#' Ysmooth <- smooth.basis(u,D,basis)
#' Y <- Ysmooth$fd
#'
#' par(mar=c(2,1,2,2),mfrow=c(1,3))
#' ftsplot(u,time,Y[1:30],space = 0.4,type=1,ylab = "",xlab = "Day",main = "Typ1=1")
#' ftsplot(u,time,Y[1:30],space = 0.4,type=2,ylab = "",xlab = "Day",main = "Typ1=2")
#' ftsplot(u,time,Y[1:30],space = 0.4,type=3,ylab = "",xlab = "Day",main = "Typ1=3")
#' @export
ftsplot = function(x, y, X, type = 2,
                   zlab = NULL, xlab = NULL, ylab = NULL,
                   space = 0.1, main = NULL) {
  u1 = as.numeric(x)
  Z = eval.fd(X, u1)
  m = nrow(Z)
  n = ncol(Z)
  zcol = matrix(rep(1:n, m),
                nrow = m, byrow = T)
  if (type == 1) {
    plot3D::ribbon3D(u1, y = y,
                     z = Z, scale = T, expand = 0.5,
                     bty = "g", phi = 25,
                     colvar = zcol, shade = 0,
                     ltheta = 320, theta = 65,
                     space = space, ticktype = "detailed",
                     d = 3, curtain = T,
                     xlab = ylab, ylab = xlab,
                     zlab = zlab, colkey = F,
                     main = main, lighting = F)
  }
  if (type == 2) {
    plot3D::ribbon3D(u1, y = y,
                     z = Z, scale = T, expand = 0.5,
                     bty = "g", phi = 15,
                     colvar = Z, shade = 0,
                     ltheta = 320, theta = 65,
                     space = space, ticktype = "detailed",
                     d = 3, curtain = F,
                     xlab = ylab, ylab = xlab,
                     zlab = zlab, colkey = F,
                     main = main, lighting = F)
  }
  if (type == 3)
    plot3D::image2D(t(Z), x = y,
                    y = x, lphi = 0, clab = zlab,
                    main = main, colkey = F,
                    xlab = xlab, ylab = ylab,
                    col = rev(grDevices::heat.colors(50)))
}
