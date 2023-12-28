#' Plot Functional Singular Spectrum Analysis Objects
#'
#' This function is a plotting method for objects of class functional singular spectrum analysis (\code{\link{fssa}}).
#' It aids users in making decisions during the grouping stage of univariate or multivariate functional singular spectrum analysis.
#'
#' @param x An object of class \code{\link{fssa}}.
#' @param d An integer representing the number of elementary components to plot.
#' @param idx A vector of indices specifying which eigen elements to plot.
#' @param idy A second vector of indices of eigen elements to plot (for \code{type="paired"}).
#' @param groups A list or vector of indices determining the grouping used for decomposition (for \code{type="wcor"}).
#' @param lwd A vector of line widths.
#' @param contrib A logical value. If \code{TRUE} (default), it displays the component's contribution to the total variance.
#' @param type The type of plot to be displayed. Possible types include:
#'   \itemize{
#'     \item \code{"values"} - Plot the square-root of eigen values (default).
#'     \item \code{"paired"} - Plot pairs of right singular function's coefficients (useful for detecting periodic components).
#'     \item \code{"wcor"} - Plot the W-correlation matrix for the reconstructed objects.
#'     \item \code{"vectors"} - Plot the right singular vectors (useful for detecting period length).
#'     \item \code{"lcurves"} - Plot left singular functions (useful for detecting period length).
#'     \item \code{"lheats"} - Heatmap plot of eigenfunctions, usable for \code{\link{funts}} variables observed over one or two-dimensional domains (useful for detecting meaningful patterns).
#'     \item \code{"periodogram"} - Periodogram plot right singular vectors (useful for detecting the frequencies of oscillations in functional data).
#'   }
#' @param vars A numeric value specifying the variable number (used in plotting MFSSA \code{"lheats"} or \code{"lcurves"}).
#' @param ylab A character vector representing the names of variables.
#' @param main The main plot title.
#' @param ... Additional arguments to be passed to methods, such as graphical parameters.
#' @seealso \code{\link{fssa}}, \code{\link{plotly_funts}}
#'
#' @examples
#' data("Callcenter")
#' L <- 28
#' U <- fssa(Callcenter, L)
#' plot(U, type = "values", d = 10)
#' plot(U, type = "vectors", d = 4)
#' plot(U, type = "paired", d = 6)
#' plot(U, type = "lcurves", d = 4, vars = 1)
#' plot(U, type = "lheats", d = 4)
#' plot(U, type = "wcor", d = 10)
#'
#'
#' @export
plot.fssa <- function(x, d = length(x$values),
                      idx = 1:d, idy = idx + 1, contrib = TRUE,
                      groups = as.list(1:d), lwd = 2,
                      type = "values", vars = NULL, ylab = NA, main = NA, ...) {
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
    if (is.na(main)) main <- "Singular Values"
    val <- sqrt(x$values)[idx]
    data_df <- data.frame(idx, val)
    p1 <- xyplot(val ~ idx,
      data = data_df, type = "o", lwd = lwd, col = "dodgerblue3", pch = 19, cex = 1.2,
      scales = list(cex.axis = 1.7, cex.main = 2, cex.lab = 1.8, cex = 0.8),
      main = main, xlab = "Components", ylab = "norms", grid = TRUE
    )
    graphics::plot(p1)
  } else if (type == "wcor") {
    W <- fwcor(x, groups)
    wplot(W, main = main)
  } else if (type == "lheats" || type == "lcurves") {
    if (is.null(vars)) vars <- 1:p
    for (j in 1:length(vars)) {
      if (type == "lheats") {
        z <- NULL
        for (i in idx) {
          if (class(x[[i]])[[1]] != "list") {
            n <- nrow(x[[i]])
            z <- c(z, as.vector(t(x[[i]])))
          } else {
            n <- nrow(x[[i]][[vars[j]]])
            z <- c(z, as.vector(t(x[[i]][[vars[j]]])))
          }
        }
        D0 <- expand.grid(
          x = 1L:L,
          y = 1L:n,
          groups = idx
        )
        D0$z <- z
        D0$groups <- factor(rep(main1,
          each = L * n
        ), levels = main1)
        if (is.na(main)) main <- "Singular functions"
        if (p > 1) {
          main <- paste(
            main, "of variable",
            ifelse(is.na(ylab), vars[j], ylab)
          )
        }
        p1 <- levelplot(
          z ~ x *
            y | groups,
          data = D0, par.strip.text = list(cex = 1.2),
          colorkey = list(labels = list(cex = 1.2)), cuts = 50L,
          xlab = "", ylab = "",
          scales = list(
            x = list(cex = c(1.2, 1.2)),
            y = list(
              cex = c(1.2, 1.2), # increase font size
              alternating = 1, # axes labels left/bottom
              tck = c(1, 0)
            )
          ), as.table = TRUE,
          main = list(main, cex = 1.5),
          col.regions = grDevices::heat.colors(100),
          ...
        )
        main <- NA
        graphics::plot(p1)
      } else if (type == "lcurves" && x$Y$dimSupp[[vars[j]]] == 1) {
        z <- NULL
        for (i in idx) {
          if (class(x[[i]])[[1]] != "list") {
            n <- nrow(x[[i]])
            z <- c(z, as.vector((x[[i]])))
          } else {
            n <- nrow(x[[i]][[vars[j]]])
            z <- c(z, as.vector((x[[i]][[vars[j]]])))
          }
        }
        if (is.na(main)) main <- "Singular functions"
        if (p > 1) {
          main <- paste(
            main, "of variable",
            ifelse(is.na(ylab), vars[j], ylab)
          )
        }
        D0 <- data.frame(z = z, x = 1:n)
        D0$curves <- rep(as.character(1:L), each = n)
        D0$groups <- factor(rep(main1, each = L * n), levels = main1)
        p1 <- xyplot(z ~ x | groups,
          group = D0$curves, lwd = lwd,
          type = "l", data = D0, par.strip.text = list(cex = 1.2),
          cuts = 50L, xlab = "", ylab = "",
          scales = list(
            x = list(cex = c(1.2, 1.2)),
            y = list(
              cex = c(1.2, 1.2), # increase font size
              alternating = 1, # axes labels left/bottom
              tck = c(1, 0)
            )
          ),
          as.table = TRUE,
          main = list(main, cex = 1.5), ...
        )

        graphics::plot(p1)
      } else if (type == "lcurves" && x$Y$dimSupp[[vars[j]]] == 2) {
        warning("The `lcurves` type for 2-dimentional fssa is not developed.")
      }
    }
  } else if (type == "vectors") {
    if (is.na(main)) main <- "Right Singular vectors"
    x0 <- c(apply(x$RVectrs[, idx], 2, scale, center = F))
    D0 <- data.frame(
      x = x0,
      time = rep(1L:K, d_idx)
    )
    D0$groups <- factor(rep(main1,
      each = K
    ), levels = main1)
    p1 <- xyplot(
      x ~ time |
        groups,
      lwd = lwd, par.strip.text = list(cex = 1.5),
      data = D0, xlab = "",
      ylab = "", main = list(main, cex = 2),
      scales = list(
        x = list(cex = c(1.4, 1.4)),
        y = list(
          cex = c(1.4, 1.4), # increase font size
          alternating = 1, # axes labels left/bottom
          tck = c(1, 0)
        )
      ),
      as.table = TRUE, type = "l",
      ...
    )
    graphics::plot(p1)
  } else if (type == "paired") {
    if (is.na(main)) main <- "Paired Plot of Right Singular Vectors"
    d_idy <- length(idy)
    if (d_idx != d_idy) stop("The length of idx and idy must be same")
    x0 <- c(apply(x$RVectrs[, idx[1]:idx[(d_idx - 1)]], 2, scale, center = F))
    y0 <- c(apply(x$RVectrs[, idy[1]:idy[(d_idy - 1)]], 2, scale, center = F))
    D0 <- data.frame(x = x0, y = y0)
    main3 <- paste(as.character(idx[1]:idx[(d_idx - 1)]), "vs", as.character(idy[1]:idy[(d_idy - 1)]))
    D0$groups <- factor(rep(main3, each = K), levels = main3)
    p1 <- xyplot(x ~ y | groups,
      data = D0, xlab = "", par.strip.text = list(cex = 1.4), lwd = lwd,
      ylab = "", main = list(main, cex = 2.0),
      scales = list(
        x = list(cex = c(1.4, 1.4)),
        y = list(
          cex = c(1.4, 1.4), # increase font size
          alternating = 1, # axes labels left/bottom
          tck = c(1, 0)
        )
      ),
      as.table = TRUE, type = "l",
      ...
    )
    graphics::plot(p1)
  } else if (type == "periodogram") {
    if (is.na(main)) main <- "Periodogram"
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
    p1 <- xyplot(
      x ~ time |
        groups,
      data = D0, xlab = "", lwd = lwd,
      ylab = "", main = main,
      scales = list(y = list(at = NULL, relation = "same")),
      as.table = TRUE, type = "l",
      ...
    )
    graphics::plot(p1)
  } else {
    stop("Unsupported type of fssa plot!!")
  }
}
