#--------------------------------------------------------------
#' Functional Time Series visualization tools using plotly
#'
#' Novel different type of plots to visualize the functional time series data objects.
#'
#' @param Y univariate or multivariate functional time series
#' @param type specifies type of plot to create
#' @param npts number of points to evaluate functional object at
#' @param main The main title.
#' @param ylab The character vector of name of variables
#' @param xlab The label of the functions arguments.
#' @param tlab The time label
#' @param var an integer Specify the variable number.
#' @importFrom plotly plot_ly add_lines layout subplot add_surface
#' @import dplyr
#' @examples
#' library(Rfssa)
#' library(fda)
#' data(Callcenter) # Read data
#' u=seq(0,1,length.out=240) # Define domain of functional data
#' d=12 # number of basis elements
#' basis=create.bspline.basis(rangeval = c(0,1),nbasis = d) # create basis object
#' Y=fts(smooth.basis(argvals = u, matrix(nrow=240,ncol=365,Callcenter$calls), basis)$fd) # create functional time series
#' plot(Y,type = "heatmap")
#' plot(Y,type = "line",var = 1)
#' plot(Y,type = "3Dsurface",var = 1)
#' plot(Y,type = "3Dline", var = 1)
#' @export
plot.fts <- function(Y,npts=100,type="line",main=NULL,ylab=NULL,xlab=NULL,tlab=NULL,var=NULL){
  p <- Y$p
  N <- Y$N
  time <- Y$time
  d <- Y$d
  u <- seq(Y$rangeval[1],Y$rangeval[2],length.out = npts)
  Pl <- list()
  if(type=="line") {
    if(is.null(var)){
      for(i in 1:p) {
        y <- as.tbl(data.frame(y=c(eval.fd(Y[[i]],u))))
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
      y <- as.tbl(data.frame(y=c(eval.fd(Y[[var]],u))))
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
        z0 <- eval.fd(Y[[i]],u)
        Pl[[i]] <- plot_ly(z = z0, x=time, y = u, type = "heatmap",
                           showscale =FALSE) %>%
          layout(yaxis = list(title = y_var),xaxis = list(title = tlab))
      }
      Pl2 <- subplot(Pl, nrows = ceiling(sqrt(p)), shareX = TRUE,
                     titleY = TRUE,titleX = TRUE) %>%
        layout(title = main)
      print(Pl2)
    } else {
      if(var > p ) var <- p
      z0 <- eval.fd(Y[[var]],u)
      if(is.null(ylab)) y_var <- paste("Variable",var) else y_var <- ylab[var]
      plot_ly(z = z0, x=time, y = u, type = "heatmap",
                         showscale =FALSE) %>%
        layout(yaxis = list(title = y_var))
    }
  } else if(type=="3Dsurface"){
    if(is.null(var) | p==1) var <- 1
    if(var>p) var <- p
    z0 <- eval.fd(Y[[var]],u)
    axx <-axy <-axz <- list(
      gridcolor="rgb(255,255,255)",
      zerolinecolor="rgb(255,255,255"
    )
    axx$title <- tlab
    axy$title <- xlab
    if(is.null(ylab)) y_var <- paste("Variable",var) else y_var <- ylab[var]
    axz$title <- y_var
    plot_ly(z = z0, x=time, y = u) %>%
      layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz))%>%
      add_surface(showscale=FALSE)
  } else if(type=="3Dline"){
    if(is.null(var) | p==1) var <- 1
    if(var>p) var <- p
    D0 <- as.tbl(data.frame(z=c(eval.fd(Y[[var]],u))))
    D0$time <- as.character(rep(time,each=npts))
    D0$x <- rep(u,length = npts)
    D0 %>%
      group_by(time) %>%
      plot_ly(x=~x,z=~z,y=~time, type = 'scatter3d', mode = 'lines',
              line = list(width = 4, color = ~z, colorscale = list(c(0,'#BA52ED'), c(1,'#FCB040')))) %>%
      layout(yaxis = list(title = paste("Variable",var)))
    } else stop("The type for the plot is not valid.")
}

