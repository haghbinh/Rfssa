#--------------------------------------------------------------
#' Functional Time Series Visualization Tools Using Plotly
#'
#' This is a plotting method for univariate or multivariate functional time series (\code{\link{fts}}). This method is designed to help the user visualize
#' \code{\link{fts}} data using a variety of techniques that use plotly.
#'
#' @param x An object of class \code{\link{fts}}.
#' @param vars A numeric specifying which variables in the fts to plot. The default is to plot all variables in succession. Note as well that variable indices may be repeated.
#' @param types A tuple of strings specifying the types of plots to be displayed where possible types for fts variables observed over a one-dimensional domain are:
#' \itemize{
#' \item \code{"line"} plot the \code{\link{fts}} elements in a line plot (default)
#' \item \code{"heatmap"} plot the \code{\link{fts}} elements in a heat map which can be used for variables observed over one or two-dimensional domains
#' \item \code{"3Dsurface"} plot the \code{\link{fts}} elements as a surface
#' \item \code{"3Dline"} plot the \code{\link{fts}} elements in a three-dimensional line plot.
#' }
#' The current plot type supported for fts variables observed over a two-dimensional domain is \code{"heatmap"}. Also note that
#' the same variable may be plotted several times using many different type options.
#' @param subplot A logical specifying whether or not line plots should be plotted in a subplot or not. The default is \code{TRUE} and if any other plot type is provided, the value is switched to \code{FALSE}.
#' @param mains A tuple of strings providing the the main titles of each plot.
#' @param ylabels A tuple of strings providing the the y-axis titles of each plot.
#' @param xlabels A tuple of strings providing the the x-axis titles of each plot.
#' @param tlabels A tuple of strings providing the the time-axis titles of each plot.
#' @param zlabels A tuple of strings providing the the z-axis titles of each plot.
#' @param ... arguments to be passed to methods, such as graphical parameters.
#' @importFrom plotly plot_ly add_lines layout subplot add_surface hide_colorbar ggplotly
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_distiller xlab ylab labs ggtitle
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom tibble as_tibble
#' @import dplyr
#'
#' @note For examples, see \code{\link{fssa}}
#'
#' @export
plot.fts <- function(x, vars = NULL, types = NULL, subplot = TRUE, mains = NULL, ylabels = NULL, xlabels = NULL, tlabels = NULL, zlabels = NULL, ...) {
  p <- length(x@C)
  N <- ncol(x@C[[1]])
  time <- colnames(x@C[[1]])
  count_twod <- 0
  Pl <- list()

  if (is.null(types)) types <- rep(NA, p)
  if (is.null(ylabels)) ylabels <- rep(NA, p)
  if (is.null(xlabels)) xlabels <- rep(NA, p)
  if (is.null(tlabels)) tlabels <- rep(NA, p)
  if (is.null(zlabels)) zlabels <- rep(NA, p)
  if (is.null(mains)) mains <- rep(NA, p)
  if (is.null(vars) == FALSE && length(types) != length(vars)) warning("\"vars\" and \"types\" are not the same length. Some plots might not appear as expected.")


  if (is.null(vars) == TRUE) vars <- 1:p
  for (j in 1:length(vars)) {
    if (ncol(x@grid[[vars[j]]]) == 1) {
      if (is.na(types[j]) == TRUE || types[j] == "line") {
        if (is.na(ylabels[j])) ylabels[j] <- "y"
        if (is.na(xlabels[j])) xlabels[j] <- "x"
        if (is.na(mains[j]) && length(vars) == 1 || is.na(mains[j]) && subplot == FALSE) mains[j] <- paste("Variable", vars[j])
        if (subplot == TRUE && length(vars) > 1) mains[j] <- NA
        y <- tibble::as_tibble(data.frame(y = c(x@B[[vars[j]]] %*% x@C[[vars[j]]])))
        y$time <- as.factor(rep(time, each = nrow(x@grid[[vars[j]]])))
        y$x <- rep(1:nrow(x@grid[[vars[j]]]), ncol(x@C[[vars[j]]]))
        Pl[[j]] <- y %>%
          group_by(time) %>%
          plot_ly(x = ~x, y = ~y) %>%
          add_lines(color = ~time, colors = c("lightsteelblue", "royalblue4"), showlegend = FALSE) %>%
          layout(title = mains[j], yaxis = list(title = ylabels[j]), xaxis = list(title = xlabels[j]))
      } else if (types[j] == "heatmap") {
        u <- seq(x@grid[[vars[j]]])
        if (is.na(xlabels[j])) xlabels[j] <- "x"
        if (is.na(tlabels[j])) tlabels[j] <- "t"
        if (is.na(mains[j]) && length(vars) == 1 || is.na(mains[j]) && subplot == FALSE) mains[j] <- paste("Variable", vars[j])
        if (subplot == TRUE && length(vars) > 1) mains[j] <- NA
        z0 <- x@B[[vars[j]]] %*% x@C[[vars[j]]]
        Pl[[j]] <- plot_ly(
          z = z0, x = time, y = u, type = "heatmap", colorscale = list(c(0, "#FFFFFAFF"), c(1, "#FF0000FF")),
          showscale = FALSE
        ) %>%
          layout(title = mains[j], yaxis = list(title = xlabels[j]), xaxis = list(title = tlabels[j]))
      } else if (types[j] == "3Dsurface") {
        u <- seq(x@grid[[vars[j]]])
        z0 <- x@B[[vars[j]]] %*% x@C[[vars[j]]]
        axx <- axy <- axz <- list(
          gridcolor = "rgb(180, 180, 180)",
          zerolinecolor = "rgb(255,255,255)"
        )
        axx$title <- ifelse(is.na(tlabels[j]), "t", tlabels[j])
        axy$title <- ifelse(is.na(xlabels[j]), "x", xlabels[j])
        axz$title <- ifelse(is.na(zlabels[j]), paste("Variable", vars[j]), zlabels[j])
        Pl[[j]] <- plot_ly(z = z0, x = time, y = u, colorscale = list(c(0, "#FFFFFAFF"), c(1, "#FF0000FF"))) %>%
          layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz)) %>%
          add_surface(showscale = FALSE)
      } else if (types[j] == "3Dline") {
        D0 <- tibble::as_tibble(data.frame(z = c(x@B[[vars[j]]] %*% x@C[[vars[j]]])))
        D0$time <- rep(time, each = nrow(x@grid[[vars[j]]]))
        D0$x <- rep(1:nrow(x@grid[[vars[j]]]), ncol(x@C[[vars[j]]]))
        axx <- axy <- axz <- list(
          gridcolor = "rgb(180, 180, 180)",
          zerolinecolor = "rgb(255,255,255)"
        )
        axx$title <- ifelse(is.na(tlabels[j]), "t", tlabels[j])
        axy$title <- ifelse(is.na(xlabels[j]), "x", xlabels[j])
        axz$title <- ifelse(is.na(zlabels[j]), paste("Variable", vars[j]), zlabels[j])
        Pl[[j]] <- D0 %>%
          group_by(time) %>%
          plot_ly(
            x = ~time, z = ~z, y = ~x, type = "scatter3d", mode = "lines", color = ~z,
            line = list(width = 4), colors = c("#FFFFFAFF", "#FF0000FF")
          ) %>%
          layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz)) %>%
          hide_colorbar()
      }
    } else {
      if (is.na(types[j]) == FALSE && types[j] != "heatmap") warning("The only plotting option available for variables observed over two-dimensional domains is \"heatmap\". Other plotting options will be added for these types of variables in the future.")
      count_twod <- count_twod + 1
      x_1 <- NULL
      x_2 <- NULL
      time <- NULL
      if (is.na(ylabels[j])) ylabels[j] <- "y"
      if (is.na(xlabels[j])) xlabels[j] <- "x"
      if (is.na(zlabels[j])) zlabels[j] <- "z"
      if (is.na(mains[j])) mains[j] <- paste("Variable", vars[j])
      y <- tibble::as_tibble(data.frame(y = c(x@B[[vars[j]]] %*% x@C[[vars[j]]])))
      time <- as.character(1:ncol(x@C[[1]]))
      y$time <- as.factor(rep(time, each = nrow(x@grid[[vars[j]]])))
      y$x_1 <- rep(x@grid[[vars[j]]][, 1], ncol(x@C[[vars[j]]]))
      y$x_2 <- rep(x@grid[[vars[j]]][, 2], ncol(x@C[[vars[j]]]))
      Pl[[j]] <- ggplotly(ggplot(y, aes(x_1, x_2, fill = y, frame = time)) +
        geom_tile() +
        scale_fill_distiller(palette = "RdYlBu") +
        theme_ipsum() +
        xlab(xlabels[j]) +
        ylab(ylabels[j]) +
        labs(fill = zlabels[j]) +
        ggtitle(mains[j]))
    }
  }

  if (count_twod >= 1) {
    print(Pl)
  } else {
    if ("3Dsurface" %in% types || "3Dline" %in% types) {
      print(Pl)
    } else if ("3Dsurface" %in% types == FALSE && subplot != TRUE || "3Dline" %in% types == FALSE && subplot != TRUE) {
      print(Pl)
    } else {
      print(subplot(Pl))
    }
  }
}
