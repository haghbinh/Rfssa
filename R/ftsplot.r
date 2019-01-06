#--------------------------------------------------------------
#' 3D fts Plots with Plot3D
#'
#' Define some new plots for functional time series data.
#' @param x a numeric vector of cordinate x.
#' @param y a numeric vector of cordinate y.
#' @param X a functional time series object of class "fd".
#' @param type what type of plot should be drawn. Possible types are
#'     "type=1" for ribbon3D plot with curtain, "type=2" for ribbon3D plot
#'      without curtain, "type=3" for image2D plot.
#' @param zlab The lable of z axis.
#' @param xlab The lable of x axis.
#' @param ylab The lable of y axis.
#' @param space The amount of space (as a fraction of the average ribbon width) left between ribbons.
#' @param main The main title.
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
                     z = Z, scale = T, expand = 0.6,
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
