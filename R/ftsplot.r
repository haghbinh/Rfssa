#--------------------------------------------------------------
#' Functional Time Series visualization tools using plotly
#'
#' Novel different type of plots to visualize the functional time series data objects.
#'
#' @param Y ??
#' @param type ??
#' @param npts ??
#' @param type ??
#' @param main The main title.
#' @param ylab The character vector of name of variables
#' @param xlab The label of the functions arguments.
#' @param tlab The time label
#' @param var an integer Specify the variable number.
#' @importFrom plotly plot_ly add_lines layout subplot add_surface
#' @importFrom fda eval.fd
#' @examples
#' plot(Y,type = "heat")
#' plot(Y,type = "3D",var = 1)
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
        y <- as.tbl(data.frame(y=c(eval.fd(Y$fd[[i]],u))))
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
      y <- as.tbl(data.frame(y=c(eval.fd(Y$fd[[var]],u))))
      y$time <- as.factor(rep(time,each=npts))
      y$x <- rep(u,length = npts)
      Pl <- y %>%
        group_by(time) %>%
        plot_ly(x=~x,y=~y) %>%
        add_lines(color = ~time,colors=c("lightsteelblue1","royalblue4"),
                  showlegend=FALSE) %>%
        layout(yaxis = list(title = paste("Variable",var)))
      print(Pl)
    }

  } else if(type=="heatmap")  {
    if(is.null(var)){
      for(i in 1:p) {
        if(is.null(ylab)) y_var <- paste("Variable",i) else y_var <- ylab[i]
        z0 <- eval.fd(Y$fd[[i]],u)
        Pl[[i]] <- plot_ly(z = z0, x=time, y = u, type = "heatmap",
                           showscale =FALSE) %>%
          layout(yaxis = list(title = y_var),xaxis = list(title = tlab))
      }
      Pl2 <- subplot(Pl, nrows = ceiling(sqrt(p)), shareX = TRUE,
                     titleY = TRUE,titleX = TRUE) %>%
        layout(title = main)
      print(Pl2)
    } else {
      z0 <- eval.fd(Y$fd[[var]],u)
      if(is.null(ylab)) y_var <- paste("Variable",var) else y_var <- ylab[var]
      Pl <- plot_ly(z = z0, x=time, y = u, type = "heatmap",
                         showscale =FALSE) %>%
        layout(yaxis = list(title = y_var))
      print(Pl)
    }
  } else if(type=="3D"){
    if(is.null(var)) var <- 1
    z0 <- eval.fd(Y$fd[[var]],u)
    axx <-axy <-axz <- list(
      gridcolor="rgb(255,255,255)",
      zerolinecolor="rgb(255,255,255"
    )
    axx$title <- tlab
    axy$title <- xlab
    if(is.null(ylab)) y_var <- paste("Variable",var) else y_var <- ylab[var]
    axz$title <- y_var
    Pl <- plot_ly(z = z0, x=time, y = u) %>%
      layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz))%>%
      add_surface(showscale=FALSE)
    print(Pl)
  } else stop("The type for the plot is not valid.")
}

