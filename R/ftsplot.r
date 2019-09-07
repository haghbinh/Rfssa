#--------------------------------------------------------------
#' Functional Time Series Visualization Tools Using Plotly
#'
#' This is a plotting method for univariate or multivariate functional time series (\code{\link{fts}}). This method is designed to help the user visualize
#' \code{\link{fts}} data using a variety of techniques that use plotly.
#'
#' @param x an object of class \code{\link{fts}}
#' @param type the type of plot to be displayed where possible types are:
#' \itemize{
#' \item \code{"line"} plot the \code{\link{fts}} elements in a line plot (default)
#' \item \code{"heatmap"} plot the \code{\link{fts}} elements in a heat map
#' \item \code{"3Dsurface"} plot the \code{\link{fts}} elements as a surface
#' \item \code{"3Dline"} plot the \code{\link{fts}} elements in a three-dimensional line plot
#' }
#' @param npts number of points to evaluate functional object at
#' @param main the main title
#' @param ylab the y-axis label
#' @param xlab the x-axis label
#' @param tlab the time-axis label
#' @param var an integer specifying the variable number to plot if \code{type="3Dsurface"} or \code{type="3Dline"}
#' @param ... arguments to be passed to methods, such as graphical parameters.
#' @importFrom plotly plot_ly add_lines layout subplot add_surface
#' @import dplyr
#' @examples
#'
#' \dontrun{
#' require(fda)
#' require(Rfssa)
#' data(Callcenter) # Read data
#' u=seq(0,1,length.out=240) # Define domain of functional data
#' d=12 # number of basis elements
#' basis=create.bspline.basis(rangeval = c(0,1),nbasis = d) # create basis object
#' smooth.calls=smooth.basis(u, matrix(nrow=240,ncol=365,Callcenter$calls), basis)
#' Y=fts(smooth.calls$fd) # create functional time series
#' plot(Y,type = "heatmap")
#' plot(Y,type = "line",var = 1)
#' plot(Y,type = "3Dsurface",var = 1)
#' }
#'
#' @note for a multivariate example, see the examples in \code{\link{fssa}}
#'
#' @references
#' Carson Sievert (2018) plotly for R. https://plotly-r.com
#'
#' @export
plot.fts <- function(x,npts=100,type="line",main=NULL,ylab=NULL,xlab=NULL,tlab=NULL,var=NULL, ...){
  p <- x$p
  N <- x$N
  time <- x$time
  d <- x$d
  u <- seq(x$rangeval[1],x$rangeval[2],length.out = npts)
  Pl <- list()
  if(type=="line") {
    if(is.null(var)){
      for(i in 1:p) {
        y <- as.tbl(data.frame(y=c(eval.fd(x[[i]],u))))
        y$time <- as.factor(rep(time,each=npts))
        y$x <- rep(u,length = npts)
        if(is.null(ylab)) y_var <- paste("Variable",i) else y_var <- ylab[i]
        Pl[[i]] <- y %>%
          group_by(time) %>%
          plot_ly(x=~x,y=~y) %>%
          add_lines(color = ~time,colors=c("lightsteelblue1","royalblue4"),
                    showlegend=FALSE) %>%
          layout(yaxis = list(title = y_var),xaxis = list(title = xlab))
      }
      Pl2 <- subplot(Pl, nrows = ceiling(sqrt(p)), shareX = TRUE,
                     titleY = TRUE,titleX = TRUE) %>%
        layout(title = main)
      print(Pl2)
    } else {
      if(var>p) var <- p
      if(is.null(ylab)) y_var <- paste("Variable",var) else y_var <- ylab[var]
      y <- as.tbl(data.frame(y=c(eval.fd(x[[var]],u))))
      y$time <- as.factor(rep(time,each=npts))
      y$x <- rep(u,length = npts)
      y %>%
        group_by(time) %>%
        plot_ly(x=~x,y=~y) %>%
        add_lines(color = ~time,colors=c("lightsteelblue1","royalblue4"),
                  showlegend=FALSE) %>%
        layout(yaxis = list(title = y_var))
    }

  } else if(type=="heatmap")  {
    if(is.null(var)){
      for(i in 1:p) {
        if(is.null(ylab)) y_var <- paste("Variable",i) else y_var <- ylab[i]
        z0 <- eval.fd(x[[i]],u)
        Pl[[i]] <- plot_ly(z = z0, x=time, y = u, type = "heatmap", colorscale = list(c(0,'#FFFFFAFF'), c(1,'#FF0000FF')),
                           showscale =FALSE) %>%
          layout(yaxis = list(title = y_var),xaxis = list(title = tlab))
      }
      Pl2 <- subplot(Pl, nrows = ceiling(sqrt(p)), shareX = TRUE,
                     titleY = TRUE,titleX = TRUE) %>%
        layout(title = main)
      print(Pl2)
    } else {
      if(var > p ) var <- p
      z0 <- eval.fd(x[[var]],u)
      if(is.null(ylab)) y_var <- paste("Variable",var) else y_var <- ylab[var]
      plot_ly(z = z0, x=time, y = u, type = "heatmap", colorscale = list(c(0,'#FFFFFAFF'), c(1,'#FF0000FF')),
                         showscale =FALSE) %>%
        layout(yaxis = list(title = y_var))
    }
  } else if(type=="3Dsurface"){
    if(is.null(var) | p==1) var <- 1
    if(var>p) var <- p
    z0 <- eval.fd(x[[var]],u)
    axx <-axy <-axz <- list(
      gridcolor="rgb(196, 196, 196)",
      zerolinecolor="rgb(255,255,255)"
    )
    axx$title <- ifelse(is.null(tlab),"time",tlab)
    axy$title <- ifelse(is.null(xlab),"x",xlab)
    axz$title <- ifelse(is.null(ylab),paste("Variable",var),ylab[var])
    plot_ly(z = z0, x=time, y = u, colorscale = list(c(0,'#FFFFFAFF'), c(1,'#FF0000FF'))) %>%
      layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz))%>%
      add_surface(showscale=FALSE)
  } else if(type=="3Dline"){
    if(is.null(var) | p==1) var <- 1
    if(var>p) var <- p
    D0 <- as.tbl(data.frame(z=c(eval.fd(x[[var]],u))))
    D0$time <- as.character(rep(time,each=npts))
    D0$x <- rep(u,length = npts)
    axx <-axy <-axz <- list(
      gridcolor="rgb(196, 196, 196)",
      zerolinecolor="rgb(255,255,255)"
    )
    axx$title <- ifelse(is.null(tlab),"x",tlab)
    axy$title <- ifelse(is.null(xlab),"time",xlab)
    axz$title <- ifelse(is.null(ylab),paste("Variable",var),ylab[var])
    D0 %>%
      group_by(time) %>%
      plot_ly(y=~time,z=~z,x=~x, type = 'scatter3d', mode = 'lines', color = ~z,
              line = list(width = 4), colors=c("#FFFFFAFF","#FF0000FF")) %>%
      layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz)) %>% layout(showlegend = FALSE)
    } else stop("The type for the plot is not valid.")
}

