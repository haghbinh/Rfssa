#' Functional SSA
#'
#' fssa is a function for decomposition stage (including embeding
#'  and  SVD step) of a functional time series.
#' @return list The outputs of following function is a list which
#' includes functional eigen vectors in the form of L-variate functional object.
#' @param Y a functional time series.
#' @param L Windows length
#' @importFrom fda fd
#' @importFrom fda inprod eval.fd smooth.basis
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

#--------------------------------------------------------------
#' Plot Plot the results of FSSA Decomposition.
#'
#'  Method dispatch of an fssa objects
#' @param U a funtional singular value decomposition object
#' @param d an integer which is the number of elementary components in the plot.
#' @param type what type of plot should be drawn. Possible types are
#' \itemize{
#' \item "values" for sqruare-root of singular values plot.
#' \item "paired" for ...
#' \item "wcor" for ...
#' \item "vectors" for ...
#' \item "meanvectors" for ...
#' \item "meanpaired" for ...
#' \item "efunctions" for ...
#' \item "efunctions2" for ...
#' }

#' @export
plot.fssa <- function(U, d = length(U$values),
    type = "values") {
    val <- sqrt(U$values)[1L:d]
    A <- val/sum(val)
    pr = round(A * 100L, 2L)
    main1 = paste0(1L:d, "(", pr,"%)")
    basis <- U[[1L]]$basis
    N <- U$N
    L <- U$L
    if (type == "values") {
        plot(val, type = "o", lwd = 2L,
            col = "dodgerblue3", pch = 19L,
            cex = 0.8, main = "Singular Values",
            ylab = " ", xlab = "Components")
    }
    if (type == "paired") {
        n0 <- nrow(U$Y$coefs) * L
        x0 <- c(sapply(U[1L:d], function(x) as.vector(t(x$coefs))))
        D0 <- data.frame(x = x0[1L:((d - 1) * n0)], y = x0[(n0 + 1):(d * n0)])
        D0$group <- as.ordered(rep(paste(main1[1:(d - 1)], "vs", main1[2:d]), each = n0))
        p1 <- lattice::xyplot(x ~ y | group,
            data = D0, xlab = "",
            ylab = "", main = "Pairs of eigenvectors",
            scales = list(x = list(at = NULL),
                y = list(at = NULL)),
            as.table = TRUE, type = "l")
        plot(p1)
    }
    if (type == "wcor") {
        W = fwcor(U, d)
        wplot(W)
    }
    if (type == "vectors") {
        n0 <- nrow(U$Y$coefs) * L
        x0 <- c(sapply(U[1L:d], function(x) as.vector(t(x$coefs))))
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
        plot(p1)
    }
    if (type == "meanvectors") {
        u <- basis$rangeval
        xindx <- seq(min(u), max(u),
            length = 100)
        x0 <- c(sapply(U[1L:d],
            function(x) colMeans(eval.fd(xindx,
                x))))
        D0 <- data.frame(x = x0,
            time = rep(1L:L, d))
        D0$group <- as.ordered(rep(main1,
            each = L))
        p1 <- lattice::xyplot(x ~ time |
            group, data = D0, xlab = "",
            ylab = "", main = "Meaned Eigenveactors",
            scales = list(x = list(at = NULL),
                y = list(at = NULL)),
            as.table = TRUE, type = "l")
        plot(p1)
    }
    if (type == "meanpaired") {
        u <- basis$rangeval
        xindx <- seq(min(u), max(u),
            length = 100L)
        x0 <- c(sapply(U[1L:d],
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
            ylab = "", main = "Meaned pairs eigenvectors",
            scales = list(x = list(at = NULL),
                y = list(at = NULL)),
            as.table = TRUE, type = "l")
        plot(p1)
    }
    if (type == "efunctions") {
        u <- basis$rangeval
        xindx <- seq(min(u), max(u),
            length = 100L)
        n <- length(xindx)
        z0 <- lapply(U[1L:d], function(x) t(eval.fd(xindx,
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
            main = "Eigenfunctions",
            col.regions = heat.colors(100))
        plot(p1)
    }
    if (type == "efunctions2") {
        col2 <- rainbow(L)
        d1 <- floor(sqrt(d))
        d2 <- ifelse(d1^2 < d,
            d1 + 1L, d1)
        par(mfrow = c(d1, d2),
            mar = c(2, 2, 3, 1))
        for (i in 1:d) plot(U[[i]],
            lty = 1, xlab = "",
            main = main1[i], ylab = "",
            lwd = 2, col = col2)
        par(mfrow = c(1, 1))
    }
}
