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
#' \item \code{"values"} - plot the square-root of singular values (default)
#' \item \code{"paired"} - plot the pairs of eigenfunction's coefficients (useful for the detection of periodic components)
#' \item \code{"wcor"} - plot the W-correlation matrix for the reconstructed objects
#' \item \code{"vectors"} - plot the eigenfunction's coefficients (useful for the detection of period length)
#' \item \code{"lcurves"} - plot of the eigenfunctions (useful for the detection of period length)
#' \item \code{"lheats"} - heatmap plot of the eigenfunctions which can be used for \code{\link{fts}} variables observed over one or two-dimensional domains (useful for the detection of meaningful patterns)
#' \item \code{"periodogram"} - periodogram plot (useful for the detecting the frequencies of oscillations in functional data).
#' }
#' @param vars A numeric specifying the variable number (can be used in plotting MFSSA \code{"lheats"} or \code{"lcurves"}).
#' @param ylab The character vector of name of variables.
#' @param main The main plot title
#' @param color_palette A string specifying the color palette that is offered by the ggplot2 package to be used when plotting left singular functions corresponding with \code{\link{fts}} variables observed over two-dimensional domains.
#' @param reverse_color_palette A boolean specifying if the color palette scale should be reversed.
#' @param ... Arguments to be passed to methods, such as graphical parameters.
#' @seealso \code{\link{fssa}}, \code{\link{plot.fts}}
#' @note See \code{\link{fssa}} examples.
#' @export
plot.fssa <- function(x, d = length(x$values),
                      idx = 1:d, idy = idx + 1, contrib = TRUE,
                      groups = as.list(1:d),
                      type = "values", vars = NULL, ylab = NA, main = NA,
                      color_palette = "RdYlBu", reverse_color_palette = FALSE, ...) {
  p <- length(x$Y$dimSupp)
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
    if(is.na(main)) main = "Singular Values";
    val <- sqrt(x$values)[idx]
    graphics::plot(idx, val,
      type = "o", lwd = 3L,
      col = "dodgerblue3", pch = 19L,
      cex.axis = 1.7, cex.main = 2, cex.lab = 1.8,
      cex = 0.8, main = main,
      ylab = "norms", xlab = "Components"
    )
  } else if (type == "wcor") {
    W <- fwcor(x, groups)
    wplot(W,main=main)
  } else if (type == "lheats" || type == "lcurves") {
    flag_1 <- 0
    flag_2 <- 0
    if (is.null(vars)) vars <- 1:p

    for (j in 1:length(vars)) {
      if (type == "lheats" && ncol(x$Y$argval[[vars[j]]]) == 1) {
        flag_1 <- 1
        n <- nrow(x$Y$argval[[vars[j]]])
        z <- c(0)
        for (i in idx) {
          if (class(x[[i]])[[1]]!="list") {
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
        if(is.na(main)) main <- "Singular functions"
        if (p > 1) {
          main <- paste(
            main, "of variable",
            ifelse(is.na(ylab), vars[j], ylab)
          )
        }
        p1 <- lattice::levelplot(z ~ x *
          y | groups,
        data = D0, par.strip.text=list(cex=1.2),
        colorkey = list(labels=list(cex=1.2)), cuts = 50L,
        xlab = "", ylab = "",
        scales = list(x = list(cex=c(1.2,1.2)),
                      y = list(cex=c(1.2, 1.2), # increase font size
                               alternating=1,   # axes labels left/bottom
                               tck = c(1,0))),
        aspect = "xy", as.table = TRUE,
        main = list(main,cex=1.5),
        col.regions = grDevices::heat.colors(100)
        )
        main = NA


        graphics::plot(p1)
      } else if (type == "lheats" && ncol(x$Y$argval[[vars[j]]]) == 2) {
        flag_2 <- 1
        time <- NULL
        for (i in idx) {
          if (class(x[[i]])[[1]]!="list") {
            temp_df = data.frame(matrix(ncol = 3, nrow = ncol(x[[1]])*nrow(x[[1]])))
            colnames(temp_df)=c("z","x","y")
            temp_df$z = c(x[[i]])
            y <- tibble::as_tibble(temp_df)
          } else {
            temp_df = data.frame(matrix(ncol = 3, nrow = ncol(x[[1]][[vars[j]]])*nrow(x[[1]][[vars[j]]])))
            colnames(temp_df)=c("z","x","y")
            temp_df$z = c(x[[i]][[vars[j]]])
            y <- tibble::as_tibble(temp_df)
          }
          temp_1 = rev(rep(x$Y$argval[[vars[j]]][, 1], L))
          temp_2 = rep(x$Y$argval[[vars[j]]][, 2], L)
          y[,2] <- temp_2
          y[,3] <- temp_1
          y$Time <- as.factor(rep(1:L, each = nrow(x$Y$argval[[vars[j]]])))
          direction <- 1
          if(isTRUE(reverse_color_palette)) direction <- -1;
         print(ggplotly(ggplot(y, aes_string("x", "y", fill = "z", frame = "Time")) +
            geom_tile() +
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_line(colour = "black")) +
              scale_fill_distiller(palette = color_palette, direction = direction) +
            ggtitle(paste("Singular function", as.character(i), "of variable", as.character(vars[j]))) +
               scale_x_continuous(name = 'x') + scale_y_continuous(name = 'y')
            ))
        }

      } else if (type == "lcurves" && ncol(x$Y$argval[[vars[j]]]) == 1) {
        flag_1 <- 1
        col2 <- grDevices::rainbow(L)
        d1 <- floor(sqrt(d_idx))
        d2 <- ceiling(d_idx / d1)
        graphics::par(
          mfrow = c(d1, d2),
          mar = c(2, 2, 3, 1), oma = c(2, 2, 7, 1), cex.main = 2
        )
        if(is.na(main)) main <- "Singular functions"
        if (p > 1) {
          main <- paste(
            main, "of variable",
            ifelse(is.na(ylab), vars[j], ylab)
          )
        }

        for (i in 1:d_idx) {
          if (class(x[[i]])[[1]]!="list") {

              graphics::matplot(x$Y$argval[[j]], x[[idx[i]]],
                                type = "l",
                                lty = 1, xlab = "", ylim = range(x[[idx[i]]]),
                                main = main1[i], ylab = "",
                                lwd = 2, col = col2, cex.main = 2.7,
                                cex.axis = 2.4
              )


            graphics::title(list(main,cex=2), outer = TRUE)
            main = NA
          } else {
            graphics::matplot(x$Y$argval[[vars[j]]], x[[idx[i]]][[vars[j]]],
              type = "l",
              lty = 1, xlab = "", ylim = range(x[[idx[i]]][[vars[j]]]),
              main = main1[i], ylab = "",
              lwd = 2, col = col2, cex.main = 2.7,
              cex.axis = 2.4
            )
            graphics::title(list(main,cex=2), outer = TRUE)
            main = NA
          }
        }

        graphics::par(mfrow = c(1, 1))
      } else if (type == "lcurves" && ncol(x$Y$argval[[vars[j]]]) == 2) {
        flag_2 <- 1
        time <- NULL
        warning("\"lcurves\" plotting option defined only for variables observed over one-dimensional domains. Plotting variable singular functions using \"lheats\".")
        for (i in idx) {
          if (class(x[[i]])[[1]]!="list") {
            temp_df = data.frame(matrix(ncol = 3, nrow = ncol(x[[1]])*nrow(x[[1]])))
            colnames(temp_df)=c("z","x","y")
            temp_df$z = c(x[[i]])
            y <- tibble::as_tibble(temp_df)
          } else {
            temp_df = data.frame(matrix(ncol = 3, nrow = ncol(x[[1]][[vars[j]]])*nrow(x[[1]][[vars[j]]])))
            colnames(temp_df)=c("z","x","y")
            temp_df$z = c(x[[i]][[vars[j]]])
            y <- tibble::as_tibble(temp_df)
          }
          temp_1 = rev(rep(x$Y$argval[[vars[j]]][, 1], L))
          temp_2 = rep(x$Y$argval[[vars[j]]][, 2], L)
          y[,2] <- temp_2
          y[,3] <- temp_1
          y$Time <- as.factor(rep(1:L, each = nrow(x$Y$argval[[vars[j]]])))
          direction <- 1
          if(isTRUE(reverse_color_palette)) direction <- -1;
         print(ggplotly(ggplot(y, aes_string("x", "y", fill = "z", frame = "Time")) +
            geom_tile() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black")) +
              scale_fill_distiller(palette = color_palette, direction = direction) +
            ggtitle(paste("Singular function", as.character(i), "of variable", as.character(vars[j]))) +
              scale_x_continuous(name = 'x') + scale_y_continuous(name = 'y')))
        }

      }
    }

    if (flag_1 == 1 && flag_2 == 1) {
      warning("The plots of singular functions for variables observed over one-dimensional domains can be found in the \"plots\" tab while plots of singular functions corresponding to variables observed over two-dimensional domains can be found in the \"viewer\" tab.")
    }
  } else if (type == "vectors") {
    if(is.na(main)) main <- "Right Singular vectors"
    x0 <- c(apply(x$RVectrs[, idx], 2, scale, center = F))
    D0 <- data.frame(
      x = x0,
      time = rep(1L:K, d_idx)
    )
    D0$groups <- factor(rep(main1,
      each = K
    ), levels = main1)
    p1 <- lattice::xyplot(x ~ time |
      groups, lwd = 2, par.strip.text=list(cex=1.5),
    data = D0, xlab = "",
    ylab = "", main = list(main,cex=2),
    scales = list(x = list(cex=c(1.4,1.4)),
                  y = list(cex=c(1.4, 1.4), # increase font size
                           alternating=1,   # axes labels left/bottom
                           tck = c(1,0))),
    as.table = TRUE, type = "l"
    )
    graphics::plot(p1)
  } else if (type == "paired") {
    if(is.na(main)) main <- "Paired Plot of Right Singular Vectors"
    d_idy <- length(idy)
    if (d_idx != d_idy) stop("The length of idx and idy must be same")
    x0 <- c(apply(x$RVectrs[, idx[1]:idx[(d_idx-1)]], 2, scale, center = F))
    y0 <- c(apply(x$RVectrs[, idy[1]:idy[(d_idy-1)]], 2, scale, center = F))
    D0 <- data.frame(x = x0, y = y0)
    main3 <- paste(as.character(idx[1]:idx[(d_idx-1)]), "vs", as.character(idy[1]:idy[(d_idy-1)]))
    D0$groups <- factor(rep(main3, each = K), levels = main3)
    p1 <- lattice::xyplot(x ~ y | groups,
      data = D0, xlab = "", par.strip.text=list(cex=1.4), lwd = 2,
      ylab = "", main = list(main,cex=2.0),
      scales = list(x = list(cex=c(1.4,1.4)),
                    y = list(cex=c(1.4, 1.4), # increase font size
                             alternating=1,   # axes labels left/bottom
                             tck = c(1,0))),
      as.table = TRUE, type = "l"
    )
    graphics::plot(p1)
  } else if (type == "periodogram") {
    if(is.na(main)) main = "Periodogram";
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
    ylab = "", main = main,
    scales = list(y = list(at = NULL, relation = "same")),
    as.table = TRUE, type = "l"
    )
    graphics::plot(p1)
  } else {
    stop("Unsupported type of fssa plot!")
  }
}
