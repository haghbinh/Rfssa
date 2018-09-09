
#' Functional SSA
#'
#' fssa is a function for decomposition stage (including embeding
#'  and  SVD step) of a functional time series.
#' @return list The outputs of following function is a list which
#' includes functional eigen vectors in the form of L-variate functional object.
#' @param Y a functional time series.
#' @param L Windows length
#' @export
fssa <- function(Y, L = floor(dim(Y$coefs)[2] / 2)) {
  N <- dim(Y$coefs)[2]
  basis <- Y$basis
  d <- basis$nbasis
  K <- N - L + 1
  B <- inprod(Y,basis)
  A <- inprod(basis, basis)
  S0 <- SS(K, L,B, d)
  H <- solve(Gram(K,L,A,d))
  Q <- eigen(H%*%S0)
  Q$vectors=Re(Q$vectors)
  out <- list(NA)
  d1 <- sum(Re(Q$values) > 0.001)
  for (i in 1:(d1))
    out[[i]] = fd(Cofmat(d, L, Q$vectors[, i]), basis)
  out$values <- Re(Q$values[1:d1])
  out$L <- L
  out$N <- N
  out$Y <- Y
  class(out) <- "fssa"
  return(out)
}




#--------------------------------------------------------------
#' Reconstruction of Functional TS
#'
#' is a function for reconstruction stage (including Grouping and
#' Hankelization steps) The output is a list of functional time series corresponds to each group.
#' "U" in the input is a fssa object. "group" is a list.
#' @param U funtional object
#' @return matrix contating Cx
#' @export

freconstruct <- function(U, group = as.list(1:10)) {
  N <- U$N
  Y <- U$Y
  d <- nrow(U[[1]]$coefs)
  L <- U$L
  K <- N - L + 1
  basis <- U[[1]]$basis
  m <- length(group)
  basis <- U[[1]]$basis
  V <- inprod(basis, basis)
  out <- list()
  for (i in 1:m) {
    Cx <- matrix(NA, nr = d, nc = N)
    g <- group[[i]]
    S <- 0
    for (j in 1:length(g)) S <- S + fproj(U,g[j], d, K, L,Y)
    S <- fH(S, d)
    Cx[, 1:L] <- S[, 1,]
    Cx[, L:N] <- S[, , L]
    out[[i]] <- fd(Cx,basis)
  }
  out$values <- sqrt(U$values)
  return(out)
}




#' Wcorrelation
#'
#'  @param  U in the input is a fssa object.
#'  @param d is the number of elementary components.
#'  @export
fwcor <- function(U, d) {
  Q <- freconstruct(U,group = as.list(1:d))
  N <- U$N
  L <- ncol(U[[1]]$coefs)
  K <- N - L + 1
  w <- 1:N
  basis <- U[[1]]$basis
  G <- inprod(basis,basis)
  L1 <- min(L, K)
  K1 <- max(K, L)
  w[L1:K1] <- L1
  w[(K1 + 1):N] <- N + 1 - ((K1 + 1):N)
  out <- matrix(1, nr = d, nc = d)
  for (i in 1:(d-1)) {
    for (j in (i+1):d) {
      out[i, j] <-
        winprod(Q[[i]],Q[[j]],w,G) / sqrt(winprod(Q[[i]], Q[[i]], w,G) * winprod(Q[[j]], Q[[j]], w,G))
    }
  }
  for(i in 2:d) for(j in 1:(i-1)) out[i,j] <- out[j,i]
  return(out)
}
#' Wplot
#'
#' @param input in the input is a ...
wplot <- function(W, main = "") {
  d <- nrow(W)
  W0 <- abs(W)
  a <- min(W0); b <- max(W0-diag(1,d)); s <- sd(W0-diag(1,d))
  diag(W0) <- min(1,b+3*s)
  xylabels <- paste0("F",1:d)
  p1<- levelplot(1-W0,xlab="",ylab="",colorkey =NULL,main=paste("W-correlation matrix"),
                 scales=list(x=list(at=1:d,lab=xylabels),
                             y=list(at=1:d,lab=xylabels)),
                 col.regions=gray(seq(0, 1, length = 100)))
  plot(p1)
}




#--------------------------------------------------------------
#' 3D fts Plots
#'
#' @export
ftsplot = function(x,
                   y,
                   X,
                   type = 2,
                   zlab = NULL,
                   xlab = NULL,
                   ylab = NULL,
                   space = 0.1,
                   main = NULL) {
  u1 = as.numeric(x)
  Z = eval.fd(X, u1)
  m = nrow(Z)
  n = ncol(Z)
  zcol = matrix(rep(1:n, m), nr = m, byrow = T)
  if (type == 1) {
    plot3D::ribbon3D(
      u1,
      y = y,
      z = Z,
      scale = T,
      expand = 0.6,
      bty = "g",
      phi = 25,
      colvar = zcol,
      shade = 0.0,
      ltheta = 320,
      theta = 65,
      space = space,
      ticktype = "detailed",
      d = 3,
      curtain = T,
      xlab = ylab,
      ylab = xlab,
      zlab = zlab,
      colkey = F,
      main = main,
      lighting = F
    )
  }
  if (type == 2) {
    plot3D::ribbon3D(
      u1,
      y = y,
      z = Z,
      scale = T,
      expand = 0.5,
      bty = "g",
      phi = 15,
      colvar = Z,
      shade = 0.0,
      ltheta = 320,
      theta = 65,
      space = space,
      ticktype = "detailed",
      d = 3,
      curtain = F,
      xlab = ylab,
      ylab = xlab,
      zlab = zlab,
      colkey = F,
      main = main,
      lighting = F
    )
  }
  if (type == 3)
    plot3D::image2D(
      t(Z),
      x = y,
      y = x,
      lphi = 0,
      clab = zlab,
      main = main
      ,colkey =F,
      xlab = xlab,
      ylab = ylab,
      col = rev(heat.colors(50))
    )
}

#--------------------------------------------------------------
#' Plots of FSSA objects
#'
#' @export
plot.fssa <-
  function(U,
           d = length(U$values),
           type = "values",
           main = NULL,
           col = "dodgerblue3") {
    val <- sqrt(U$values)[1:d]
    A <- val / sum(val)
    pr = round(A * 100, 2)
    main1 = paste0(1:d, "(", pr, "%)")
    basis <- U[[1]]$basis
    N <- U$N
    L <- U$L
    if (type == "values") {
      plot(
        val,
        type = "o",
        lwd = 2,
        col = col,
        pch = 19,
        cex = 0.8,
        main = "Singular Values",
        ylab = " ",
        xlab = "Components"
      )
    }
    if (type == "paired") {
      n0 <- nrow(U$Y$coefs)*L
      x0 <- c(sapply(U[1:d],function(x) as.vector(t(x$coefs)) ))
      D0 <- data.frame(x=x0[1:((d-1)*n0)],y=x0[(n0+1):(d*n0)])
      D0 $ group <- as.ordered(rep(paste(main1[1:(d-1)],"vs",main1[2:d]),each=n0))
      p1 <- xyplot(x~y|group,data = D0,xlab="",ylab="",main="Pairs of eigenvectors",
             scales=list(x=list(at=NULL),y=list(at=NULL)),
             as.table=T,type="l")
      plot(p1)
    }
    if (type == "wcor") {
      W = fwcor(U, d)
      wplot(W, main = main)
    }
    if (type == "vectors") {
      n0 <- nrow(U$Y$coefs)*L
      x0 <- c(sapply(U[1:d],function(x) as.vector(t(x$coefs)) ))
      D0 <- data.frame(x=x0,time=rep(1:n0,d))
      D0 $ group <- as.ordered(rep(main1,each=n0))
      p1 <- xyplot(x~time|group,data = D0,xlab="",ylab="",main="Eigenvectors",
                      scales=list(x=list(at=NULL),y=list(at=NULL)),
                      as.table=T,type="l")
      plot(p1)
    }
    if (type == "meanvectors") {
      u <- basis$rangeval
      xindx <- seq(min(u), max(u), length = 100)
      x0 <- c(sapply(U[1:d],function(x) colMeans(eval.fd(xindx,x)) ))
      D0 <- data.frame(x=x0,time=rep(1:L,d))
      D0 $ group <- as.ordered(rep(main1,each=L))
      p1 <- xyplot(x~time|group,data = D0,xlab="",ylab="",main="Meaned Eigenveactors",
                   scales=list(x=list(at=NULL),y=list(at=NULL)),
                   as.table=T,type="l")
      plot(p1)
    }
    if (type == "meanpaired") {
      u <- basis$rangeval
      xindx <- seq(min(u), max(u), length = 100)
      x0 <- c(sapply(U[1:d],function(x) colMeans(eval.fd(xindx,x)) ))
      D0 <- data.frame(x=x0[1:((d-1)*L)],y=x0[(L+1):(d*L)])
      D0 $ group <- as.ordered(rep(paste(main1[1:(d-1)],"vs",main1[2:d]),each=L))
      p1 <- xyplot(x~y|group,data = D0,xlab="",ylab="",main="Meaned pairs eigenvectors",
                   scales=list(x=list(at=NULL),y=list(at=NULL)),
                   as.table=T,type="l")
      plot(p1)
    }
    if (type == "efunctions") {
      u <- basis$rangeval
      xindx <- seq(min(u), max(u), length = 100)
      n <- length(xindx)
      z0 <- lapply(U[1:d],function(x) t(eval.fd(xindx,x)) )
      z <- c(sapply(z0, function(x) as.vector(x)))
      D0 <- expand.grid(x = 1:L, y = 1:n, group = 1:d)
      D0$z <- z
      D0 $ group <- as.ordered(rep(main1,each=L*n))
      p1 <- levelplot(z~x*y|group,data=D0,colorkey =T,cuts =50,xlab="",ylab="",
                scales=list(x=list(at=NULL),y=list(at=NULL)),aspect = "xy",
                ,as.table=T,main="Eigenfunctions"
                ,col.regions=heat.colors(100))
      plot(p1)
    }
    if (type == "functions") {
      u <- basis$rangeval
      xindx = seq(min(u), max(u), length = 100)
      y = 1:L
      d1 <- floor(sqrt(d))
      d2 <- ifelse(d1 ^ 2 < d, d1 + 1, d1)
      par(mfrow = c(d1, d2), mar = c(2, 2, 3, 1))
      for (i in 1:d) {
        ftsplot(
          xindx,
          1:L,
          U[[i]],
          type = 3,
          xlab = "Time",
          ylab = "",
          main = main1[i]
        )
      }
      par(mfrow=c(1,1))
    }
    if (type == "efunctions2") {
      col2 <- rainbow(L)
      d1 <- floor(sqrt(d))
      d2 <- ifelse(d1 ^ 2 < d, d1 + 1, d1)
      par(mfrow = c(d1, d2), mar = c(2, 2, 3, 1))
      for (i in 1:d) plot(U[[i]],lty=1,xlab="",main=main1[i],ylab="",lwd=2,col=col2)
      par(mfrow=c(1,1))
    }
  }
