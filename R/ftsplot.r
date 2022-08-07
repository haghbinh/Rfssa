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
#' \item \code{"line"} - plot the \code{\link{fts}} elements in a line plot (default)
#' \item \code{"heatmap"} - plot the \code{\link{fts}} elements in a heat map which can be used for variables observed over one or two-dimensional domains
#' \item \code{"3Dsurface"} - plot the \code{\link{fts}} elements as a surface
#' \item \code{"3Dline"} - plot the \code{\link{fts}} elements in a three-dimensional line plot.
#' }
#' The current plot type supported for fts variables observed over a two-dimensional domain is \code{"heatmap"}. Also note that
#' the same variable may be plotted several times using many different type options.
#' @param subplot A logical specifying whether or not line plots should be plotted in a subplot or not. The default is \code{TRUE} and if any other plot type is provided, the value is switched to \code{FALSE}.
#' @param mains A tuple of strings providing the the main titles of each plot.
#' @param ylabels A tuple of strings providing the the y-axis titles of each plot.
#' @param xlabels A tuple of strings providing the the x-axis titles of each plot.
#' @param tlabels A tuple of strings providing the the time-axis titles of each plot.
#' @param zlabels A tuple of strings providing the the z-axis titles of each plot.
#' @param xticklabels A list of character vectors where each entry specifies the tick labels for the domain of the functions.
#' @param xticklocs A list of numerics where each entry specifies the position of the tick labels for the domain of the functions.
#' @param yticklabels A list of character vectors where each entry specifies the tick labels for the domain of the functions.
#' @param yticklocs A list of numerics where each entry specifies the position of the tick labels for the domain of the functions.
#' @param color_palette A string specifying the color palette that is offered by the ggplot2 package to be used when plotting \code{\link{fts}} variables observed over two-dimensional domains.
#' @param reverse_color_palette A boolean specifying if the color palette scale should be reversed.
#' @param ... arguments to be passed to methods, such as graphical parameters.
#' @importFrom plotly plot_ly add_lines layout subplot add_surface hide_colorbar ggplotly
#' @importFrom ggplot2 ggplot aes_string unit geom_tile scale_fill_distiller xlab ylab labs ggtitle scale_y_continuous scale_x_continuous waiver theme element_line element_text element_blank
#' @importFrom tibble as_tibble
#' @import dplyr
#'
#' @note For examples, see \code{\link{fssa}}
#'
#' @export
plot.fts <- function(x, vars = NULL, types = NULL, subplot = TRUE, mains = NULL, ylabels = NULL, xlabels = NULL, tlabels = NULL,
                     zlabels = NULL, xticklabels = NULL, xticklocs = NULL, yticklabels = NULL, yticklocs = NULL,
                     color_palette = "RdYlBu", reverse_color_palette = FALSE, ...) {
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
  if (is.null(xticklabels)) xticklabels <- as.list(rep(NA, p))
  if (is.null(xticklocs)) xticklocs <- as.list(rep(NA, p))
  if (is.null(yticklabels)) yticklabels <- as.list(rep(NA, p))
  if (is.null(yticklocs)) yticklocs <- as.list(rep(NA, p))
  if (is.null(zlabels)) zlabels <- rep(NA, p)
  if (is.null(zlabels)) zlabels <- rep(NA, p)
  if (is.null(mains)) mains <- rep(NA, p)
  if (is.null(vars) == FALSE && length(types) != length(vars)) warning("\"vars\" and \"types\" are not the same length. Some plots might not appear as expected.")

  if (is.null(vars) == TRUE) vars <- 1:p
  for (j in 1:length(vars)) {
    if (j > length(xticklocs) || j > length(xticklabels)) xticklocs[[j]] = xticklabels[[j]] = NA
    if (j > length(yticklocs) || j > length(yticklabels)) yticklocs[[j]] = yticklabels[[j]] = NA
    if ((is.na(xticklocs[[j]][1]) && is.na(xticklabels[[j]][1])==FALSE) || (is.na(xticklocs[[j]][1])==FALSE && is.na(xticklabels[[j]][1]))) warning(paste0("Please provide horizontal axis labels and label locations for variable ",as.character(j),"."))
    if ((is.na(yticklocs[[j]][1]) && is.na(yticklabels[[j]][1])==FALSE) || (is.na(yticklocs[[j]][1])==FALSE && is.na(yticklabels[[j]][1]))) warning(paste0("Please provide vertical axis labels and label locations for variable ",as.character(j),"."))
    if (ncol(x@grid[[vars[j]]]) == 1) {
      if (is.na(types[j]) == TRUE || types[j] == "line") {
        if (is.na(ylabels[j])) ylabels[j] <- "y"
        if (is.na(xlabels[j])) xlabels[j] <- "x"
        if (is.na(mains[j])) mains[j] <- paste("Variable", vars[j])
        y <- tibble::as_tibble(data.frame(y = c(x@B[[vars[j]]] %*% x@C[[vars[j]]])))
        y$time <- as.factor(rep(time, each = nrow(x@grid[[vars[j]]])))
        y$x <- rep(1:nrow(x@grid[[vars[j]]]), ncol(x@C[[vars[j]]]))
        if((is.na(xticklabels[[j]][1]) == FALSE && is.na(xticklocs[[j]][1]) == FALSE) && (is.na(yticklabels[[j]][1]) == FALSE && is.na(yticklocs[[j]][1]) == FALSE)){
          Pl[[j]] <- y %>%
            group_by(time) %>%
            plot_ly(x = ~x, y = ~y) %>%
            add_lines(color = ~time, colors = c("lightsteelblue", "royalblue4"), showlegend = FALSE) %>%
            layout(title = paste('<b>',mains[j],'</b>'),font = list(size=12), yaxis = list(title = ylabels[j],ticktext = yticklabels[[j]], tickvals = yticklocs[[j]],
                                                                                           range = list(yticklocs[[j]][1],yticklocs[[j]][length(yticklocs[[j]])])),
                   xaxis = list(title = xlabels[j],ticktext = xticklabels[[j]], tickvals = xticklocs[[j]],range = list(xticklocs[[j]][1],xticklocs[[j]][length(xticklocs[[j]])])))
        }else if((is.na(xticklabels[[j]][1]) == FALSE && is.na(xticklocs[[j]][1]) == FALSE) && (is.na(yticklabels[[j]][1])||is.na(yticklocs[[j]][1]))){
          Pl[[j]] <- y %>%
            group_by(time) %>%
            plot_ly(x = ~x, y = ~y) %>%
            add_lines(color = ~time, colors = c("lightsteelblue", "royalblue4"), showlegend = FALSE) %>%
            layout(title = paste('<b>',mains[j],'</b>'),font = list(size=12), yaxis = list(title = ylabels[j]),
                   xaxis = list(title = xlabels[j],ticktext = xticklabels[[j]], tickvals = xticklocs[[j]],range = list(xticklocs[[j]][1],xticklocs[[j]][length(xticklocs[[j]])])))
        }else if((is.na(xticklabels[[j]][1])||is.na(xticklocs[[j]][1])) && (is.na(yticklabels[[j]][1]) == FALSE && is.na(yticklocs[[j]][1]) == FALSE)){
          Pl[[j]] <- y %>%
            group_by(time) %>%
            plot_ly(x = ~x, y = ~y) %>%
            add_lines(color = ~time, colors = c("lightsteelblue", "royalblue4"), showlegend = FALSE) %>%
            layout(title = paste('<b>',mains[j],'</b>'),font = list(size=12), yaxis = list(title = ylabels[j],ticktext = yticklabels[[j]], tickvals = yticklocs[[j]],
                                                                                           range = list(yticklocs[[j]][1],yticklocs[[j]][length(yticklocs[[j]])])),
                   xaxis = list(title = xlabels[j]))
        }else{
          Pl[[j]] <- y %>%
            group_by(time) %>%
            plot_ly(x = ~x, y = ~y) %>%
            add_lines(color = ~time, colors = c("lightsteelblue", "royalblue4"), showlegend = FALSE) %>%
            layout(title = paste('<b>',mains[j],'</b>'),font = list(size=12), yaxis = list(title = ylabels[j]),
                   xaxis = list(title = xlabels[j]))
        }
      } else if (types[j] == "heatmap") {
        u <- seq(x@grid[[vars[j]]])
        if (is.na(xlabels[j])) xlabels[j] <- "x"
        if (is.na(tlabels[j])) tlabels[j] <- "t"
        if (is.na(zlabels[j])) zlabels[j] <- "z"
        if (is.na(mains[j])) mains[j] <- paste("Variable", vars[j])
        z0 <- x@B[[vars[j]]] %*% x@C[[vars[j]]]
        if(is.na(xticklabels[[j]][1])||is.na(xticklocs[[j]][1])){
          Pl[[j]] <- plot_ly(
            z = z0, x = time, y = u, type = "heatmap", colorscale = list(c(0, "#FFFFFAFF"), c(1, "#FF0000FF")),
            showscale = FALSE, hovertemplate = paste0(zlabels[j],":"," %{z}", "\n",tlabels[j],":"," %{x}", "\n",xlabels[j],":"," %{y}")
          ) %>%
            layout(title = paste('<b>',mains[j],'</b>'), font = list(size=12),yaxis = list(title = xlabels[j]),
                   xaxis = list(title = tlabels[j]))
        }else{
          Pl[[j]] <- plot_ly(
            z = z0, x = time, y = u, type = "heatmap", colorscale = list(c(0, "#FFFFFAFF"), c(1, "#FF0000FF")),
            showscale = FALSE, hovertemplate = paste0(zlabels[j],":"," %{z}", "\n",tlabels[j],":"," %{x}", "\n",xlabels[j],":"," %{y}")
          ) %>%
            layout(title = paste('<b>',mains[j],'</b>'),font = list(size=12), yaxis = list(title = xlabels[j],ticktext = xticklabels[[j]], tickvals = xticklocs[[j]]
                                                  ), xaxis = list(title = tlabels[j]))
        }
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
        if(is.na(xticklabels[[j]][1])==FALSE||is.na(xticklocs[[j]][1])==FALSE){
          axy$ticktext = xticklabels[[j]]
          axy$tickvals = xticklocs[[j]]
        }

        Pl[[j]] <- plot_ly(z = z0, x = time, y = u, colorscale = list(c(0, "#FFFFFAFF"), c(1, "#FF0000FF")),
                           hovertemplate = paste0(axz$title[j],":"," %{z}", "\n",axx$title[j],":"," %{x}", "\n", axy$title[j],":"," %{y}")) %>%
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
        if(is.na(xticklabels[[j]][1])==FALSE||is.na(xticklocs[[j]][1])==FALSE){
          axy$ticktext = xticklabels[[j]]
          axy$tickvals = xticklocs[[j]]
        }

        Pl[[j]] <- D0 %>%
          group_by(time) %>%
          plot_ly(
            x = ~time, z = ~z, y = ~x, type = "scatter3d", mode = "lines", color = ~z,
            line = list(width = 4), colors = c("#FFFFFAFF", "#FF0000FF"),hovertemplate = paste0(axz$title[j],":"," %{z}", "\n",axx$title[j],":"," %{x}", "\n", axy$title[j],":"," %{y}")
          ) %>%
          layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz)) %>%
          hide_colorbar()
      }
    } else {
      if (is.na(types[j]) == FALSE && types[j] != "heatmap") warning("The only plotting option available for variables observed over two-dimensional domains is \"heatmap\". Other plotting options will be added for these types of variables in the future.")
      count_twod <- count_twod + 1
      if (is.na(ylabels[j])) ylabels[j] <- "y"
      if (is.na(xlabels[j])) xlabels[j] <- "x"
      if (is.na(zlabels[j])) zlabels[j] <- "z"
      if (is.na(mains[j])) mains[j] <- paste("Variable", vars[j])
      if (is.na(xticklabels[[j]][[1]])) xticklabels[[j]] = waiver()
      if (is.na(xticklocs[[j]][[1]])) xticklocs[[j]] = waiver()
      if (is.na(yticklabels[[j]][[1]])) yticklabels[[j]] = waiver()
      if (is.na(yticklocs[[j]][[1]])) yticklocs[[j]] = waiver()
      temp_df = data.frame(matrix(ncol = 3, nrow = ncol(x@C[[vars[j]]])*nrow(x@grid[[vars[j]]])))
      colnames(temp_df)=c(zlabels[j],xlabels[j],ylabels[j])
      temp_df[,1] = c(x@B[[vars[j]]] %*% x@C[[vars[j]]])
      temp_df[,2] <- rep((x@grid[[vars[j]]][, 2]), ncol(x@C[[vars[j]]]))
      temp_df[,3] <- rev(rep((x@grid[[vars[j]]][, 1]), ncol(x@C[[vars[j]]])))
      y <- tibble::as_tibble(temp_df)
      # The next line is to help alleviate a bug in ggplot2
      if(time[1] %in% as.character(1:ncol(x@C[[vars[j]]]))) time <- as.numeric(time)
      y$Time <- rep(time, each = nrow(x@grid[[vars[j]]]))
      direction <- 1
      if(isTRUE(reverse_color_palette)) direction <- -1;
      Pl[[j]] <- ggplotly(ggplot(y, aes_string(xlabels[j], ylabels[j], fill = zlabels[j], frame = "Time")) +
                            geom_tile() +
                            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                  panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                            scale_fill_distiller(palette = color_palette, direction = direction) +
                            labs(fill = zlabels[j]) +
                            ggtitle(mains[j]) +
                            scale_y_continuous(name = ylabels[j],breaks = yticklocs[[j]]
                                               ,labels = yticklabels[[j]]) +
                            scale_x_continuous(name = xlabels[j],breaks = xticklocs[[j]]
                          ,labels = xticklabels[[j]])+theme(text = element_text(size=10))
                          + theme(axis.text.x=element_text(size=10),
                                 axis.title.x=element_text(size=15))
                          + theme(axis.text.y=element_text(angle=90,size=10),
                                 axis.title.y=element_text(angle=90,size=15))
                          + theme(plot.title = element_text(size = 30, face = "bold"))
                          + theme(legend.position = "right", legend.key.height = unit(0.67, "cm"), legend.title = element_text(size=15),
                                  legend.text = element_text(size=10))
                          )

    }
  }


  if (p == 1 || "heatmap" %in% types || "3Dsurface" %in% types || "3Dline" %in% types || count_twod >=1 || subplot == FALSE || length(vars) == 1) {
    print(Pl)
  } else {
    print(subplot(Pl, titleX = TRUE, titleY = TRUE)%>%
            layout(title = ""))
  }


}
