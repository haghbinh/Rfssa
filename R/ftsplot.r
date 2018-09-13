#--------------------------------------------------------------
#3D fts Plots
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
                    col = rev(heat.colors(50)))
}
