#--------------------------------------------------------------
#' Plot Functional Singular Spectrum Analysis Objects
#'
#'  This is a plotting method for objects of class functional singular spectrum analysis (\code{\link{fssa}}). The method is designed to help the user make decisions
#'  on how to do the grouping stage of univariate or multivariate functional singular spectrum analysis.
#'
#' @param x An object of class \code{\link{fssa}}.
#' @param d An integer which is the number of elementary components in the plot.
#' @param idx A vector of indices of eigen elements to plot.
#' @param idy A second vector of indices of eigen elements to plot (for \code{type="paired"}).
#' @param groups A list or vector of indices determines grouping used for the decomposition(for \code{type="wcor"}).
#' @param contrib A logical where if the value is \code{TRUE} (the default), the contribution of the component to the total variance is displayed.
#' @param type The type of plot to be displayed where possible types are:
#' \itemize{
#' \item \code{"values"} plot the square-root of singular values (default)
#' \item \code{"paired"} plot the pairs of eigenfunction's coefficients (useful for the detection of periodic components)
#' \item \code{"wcor"} plot the W-correlation matrix for the reconstructed objects
#' \item \code{"vectors"} plot the eigenfunction's coefficients (useful for the detection of period length)
#' \item \code{"lcurves"} plot of the eigenfunctions (useful for the detection of period length)
#' \item \code{"lheats"} heatmap plot of the eigenfunctions which can be used for \code{\link{fts}} variables observed over one or two-dimensional domains (useful for the detection of meaningful patterns)
#' \item \code{"periodogram"} periodogram plot (useful for the detecting the frequencies of oscillations in functional data).
#' }
#' @param vars A numeric specifying the variable number (can be used in plotting MFSSA \code{"lheats"} or \code{"lcurves"}).
#' @param ylab The character vector of name of variables.
#' @param ... Arguments to be passed to methods, such as graphical parameters.
#' @seealso \code{\link{fssa}}, \code{\link{plot.fts}}
#' @note See \code{\link{fssa}} examples.
#' @export
plot.fssa <- function(x, d = length(x$values),
                      idx = 1:d, idy = idx + 1, contrib = TRUE,
                      groups = as.list(1:d),
                      type = "values", vars = NULL, ylab = NA, ...) {
  p <- length(x$Y@C)
  A <- ((x$values) / sum(x$values))
  pr <- round(A * 100L, 2L)
  idx <- sort(idx)
  idy <- sort(idy)
  if (max(idx) > d | min(idx) < 1) stop("The idx must be subset of 1:d.")
  d_idx <- length(idx)
  if (contrib) {
    main1 <- paste0(idx, "(", pr[idx], "%)")
    main2 <- paste0(idy, "(", pr[idy], "%)")
  } else {
    main1 <- paste(idx)
    main2 <- paste(idy)
  }
  N <- x$N
  L <- x$L
  K <- N - L + 1L

  if (type == "values") {
    val <- sqrt(x$values)[idx]
    graphics::plot(idx, val,
      type = "o", lwd = 2L,
      col = "dodgerblue3", pch = 19L,
      cex = 0.8, main = "Singular Values",
      ylab = "norms", xlab = "Components"
    )
  } else if (type == "wcor") {
    W <- fwcor(x, groups)
    wplot(W)
  } else if (type == "lheats" || type == "lcurves") {
    flag_1 <- 0
    flag_2 <- 0
    if (is.null(vars)) vars <- 1:p

    for (j in 1:length(vars)) {
      if (type == "lheats" && ncol(x$Y@grid[[vars[j]]]) == 1) {
        flag_1 <- 1
        n <- nrow(x$Y@grid[[vars[j]]])
        z <- c(0)
        for (i in idx) {
          if (p == 1) {
            z <- c(z, as.vector(t(x[[i]])))
          } else {
            z <- c(z, as.vector(t(x[[i]][[vars[j]]])))
          }
        }
        z <- z[2:length(z)]
        D0 <- expand.grid(
          x = 1L:L,
          y = 1L:n, groups = idx
        )
        D0$z <- z
        D0$groups <- factor(rep(main1,
          each = L * n
        ), levels = main1)
        title0 <- "Singular functions"
        if (p > 1) {
          title0 <- paste(
            title0, "of the variable",
            ifelse(is.na(ylab), vars[j], ylab)
          )
        }
        p1 <- lattice::levelplot(z ~ x *
          y | groups,
        data = D0,
        colorkey = TRUE, cuts = 50L,
        xlab = "", ylab = "",
        scales = list(
          x = list(at = NULL),
          y = list(at = NULL)
        ),
        aspect = "xy", as.table = TRUE,
        main = title0,
        col.regions = grDevices::heat.colors(100)
        )


        graphics::plot(p1)
      } else if (type == "lheats" && ncol(x$Y@grid[[vars[j]]]) == 2) {
        flag_2 <- 1
        x_1 <- NULL
        x_2 <- NULL
        time <- NULL
        for (i in idx) {
          if (p == 1) {
            y <- tibble::as_tibble(data.frame(y = c(x[[i]])))
          } else {
            y <- tibble::as_tibble(data.frame(y = c(x[[i]][[vars[j]]])))
          }
          y$time <- as.factor(rep(1:L, each = nrow(x$Y@grid[[vars[j]]])))
          y$x_1 <- rep(x$Y@grid[[vars[j]]][, 1], L)
          y$x_2 <- rep(x$Y@grid[[vars[j]]][, 2], L)
          print(ggplotly(ggplot(y, aes(x_1, x_2, fill = y, frame = time)) +
            geom_tile() +
            scale_fill_distiller(palette = "RdYlBu") +
            theme_ipsum() +
            ggtitle(paste("Singular function", as.character(i), "of variable", as.character(vars[j])))))
        }
      } else if (type == "lcurves" && ncol(x$Y@grid[[vars[j]]]) == 1) {
        flag_1 <- 1
        col2 <- grDevices::rainbow(L)
        d1 <- floor(sqrt(d_idx))
        d2 <- ceiling(d_idx / d1)
        graphics::par(
          mfrow = c(d1, d2),
          mar = c(2, 2, 3, 1), oma = c(2, 2, 7, 1), cex.main = 1.6
        )
        title0 <- "Singular functions"
        if (p > 1) {
          title0 <- paste(
            title0, "of the variable",
            ifelse(is.na(ylab), vars[j], ylab)
          )
        }

        for (i in 1:d_idx) {
          if (p == 1) {
            graphics::matplot(x$Y@grid[[j]], x[[idx[i]]],
              type = "l",
              lty = 1, xlab = "", ylim = range(x[[idx[i]]]),
              main = main1[i], ylab = "",
              lwd = 2, col = col2
            )
            graphics::title(title0, outer = TRUE)
          } else {
            graphics::matplot(x$Y@grid[[j]], x[[idx[i]]][[vars[j]]],
              type = "l",
              lty = 1, xlab = "", ylim = range(x[[idx[i]]][[vars[j]]]),
              main = main1[i], ylab = "",
              lwd = 2, col = col2
            )
            graphics::title(title0, outer = TRUE)
          }
        }
        graphics::par(mfrow = c(1, 1))
      } else if (type == "lcurves" && ncol(x$Y@grid[[vars[j]]]) == 2) {
        flag_2 <- 1
        x_1 <- NULL
        x_2 <- NULL
        time <- NULL
        warning("\"lcurves\" plotting option defined only for variables observed over one-dimensional domains. Plotting variable singular functions using \"lheats\".")
        for (i in idx) {
          if (p == 1) {
            y <- tibble::as_tibble(data.frame(y = c(x[[i]])))
          } else {
            y <- tibble::as_tibble(data.frame(y = c(x[[i]][[vars[j]]])))
          }
          y$time <- as.factor(rep(1:L, each = nrow(x$Y@grid[[vars[j]]])))
          y$x_1 <- rep(x$Y@grid[[vars[j]]][, 1], L)
          y$x_2 <- rep(x$Y@grid[[vars[j]]][, 2], L)
          print(ggplotly(ggplot(y, aes(x_1, x_2, fill = y, frame = time)) +
            geom_tile() +
            scale_fill_distiller(palette = "RdYlBu") +
            theme_ipsum() +
            ggtitle(paste("Singular function", as.character(i), "of variable", as.character(vars[j])))))
        }
      }
    }

    if (flag_1 == 1 && flag_2 == 1) {
      warning("The plots of singular functions for variables observed over one-dimensional domains can be found in the \"plots\" tab while plots of singular functions corresponding to variables observed over two-dimensional domains can be found in the \"viewer\" tab.")
    }
  } else if (type == "vectors") {
    x0 <- c(apply(x$RVectrs[, idx], 2, scale, center = F))
    D0 <- data.frame(
      x = x0,
      time = rep(1L:K, d_idx)
    )
    D0$groups <- factor(rep(main1,
      each = K
    ), levels = main1)
    p1 <- lattice::xyplot(x ~ time |
      groups,
    data = D0, xlab = "",
    ylab = "", main = "Singular vectors",
    scales = list(
      x = list(at = NULL),
      y = list(at = NULL, relation = "same")
    ),
    as.table = TRUE, type = "l"
    )
    graphics::plot(p1)
  } else if (type == "paired") {
    d_idy <- length(idy)
    if (d_idx != d_idy) stop("The length of idx and idy must be same")
    x0 <- c(apply(x$RVectrs[, idx], 2, scale, center = F))
    y0 <- c(apply(x$RVectrs[, idy], 2, scale, center = F))
    D0 <- data.frame(x = x0, y = y0)
    main3 <- paste(main1, "vs", main2)
    D0$groups <- factor(rep(main3, each = K), levels = main3)
    p1 <- lattice::xyplot(x ~ y | groups,
      data = D0, xlab = "",
      ylab = "", main = "Paired Singular vectors (Right)",
      scales = list(
        x = list(at = NULL, relation = "same"),
        y = list(at = NULL, relation = "same")
      ),
      as.table = TRUE, type = "l"
    )
    graphics::plot(p1)
  } else if (type == "periodogram") {
    ff <- function(x) {
      I <- abs(fft(x) / sqrt(K))^2
      P <- (4 / K) * I
      return(P[1:(floor(K / 2) + 1)])
    }
    x0 <- c(apply(apply(x$RVectrs[, idx], 2, scale, center = F), 2, ff))
    D0 <- data.frame(
      x = x0,
      time = rep((0:floor(K / 2)) / K, d_idx)
    )
    D0$groups <- factor(rep(main1,
      each = (floor(K / 2) + 1)
    ), levels = main1)
    p1 <- lattice::xyplot(x ~ time |
      groups,
    data = D0, xlab = "",
    ylab = "", main = "Periodogram of Singular vectors",
    scales = list(y = list(at = NULL, relation = "same")),
    as.table = TRUE, type = "l"
    )
    graphics::plot(p1)
  } else {
    stop("Unsupported type of fssa plot!")
  }
}
