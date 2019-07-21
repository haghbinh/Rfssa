#--------------------------------------------------------------
#' Functional Time Series Plots
#'
#' Novel different type of plots to visualize the functional time series data objects.
#'
#' @param Y ??
#' @param type ??
#' @param npts ??
#' @param type ??
#' @param main The main title.
#' @param var an integer Specify the variable number.
#' @examples
#' plot(Y,type = "heat")
#' plot(Y,type = "3D",var = 1)
#' @export
plot.fts <- function(Y,npts=100,type="line",main="",var=NA){
  p <- Y$p
  N <- Y$N
  time <- Y$time
  d <- Y$d
  u <- seq(Y$rangeval[1],Y$rangeval[2],length.out = npts)
  Pl <- list()
  if(type=="line") {
    if(is.na(var)){
      for(i in 1:p) {
        y <- as.tbl(data.frame(y=c(eval.fd(Y[[i]],u))))
        y$time <- as.factor(rep(time,each=npts))
        y$x <- rep(u,length = npts)
        Pl[[i]] <- y %>%
          group_by(time) %>%
          plot_ly(x=~x,y=~y) %>%
          add_lines(color = ~time,colors=c("lightsteelblue1","royalblue4"),
                    showlegend=FALSE) %>%
          layout(yaxis = list(title = paste("Variable",i)))
      }
      Pl2 <- subplot(Pl, nrows = ceiling(sqrt(p)), shareX = TRUE,
                     titleY = TRUE) %>%
        layout(title = "main")
      print(Pl2)
    } else {
      y <- as.tbl(data.frame(y=c(eval.fd(Y[[var]],u))))
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

  } else if(type=="heat")  {
    if(is.na(var)){
      for(i in 1:p) {
        z0 <- eval.fd(Y[[i]],u)
        Pl[[i]] <- plot_ly(z = z0, x=time, y = u, type = "heatmap",
                           showscale =FALSE) %>%
          layout(yaxis = list(title = paste("Variable",i)))
      }
      Pl2 <- subplot(Pl, nrows = ceiling(sqrt(p)), shareX = TRUE,
                     titleY = TRUE) %>%
        layout(title = "main")
      print(Pl2)
    } else {
      z0 <- eval.fd(Y[[var]],u)
      Pl <- plot_ly(z = z0, x=time, y = u, type = "heatmap",
                         showscale =FALSE) %>%
        layout(yaxis = list(title = paste("Variable",var)))
      print(Pl)
    }
  } else if(type=="3D"){
    if(is.na(var)) var <- 1
    z0 <- eval.fd(Y[[var]],u)
    axx <-axy <-axz <- list(
      gridcolor="rgb(255,255,255)",
      zerolinecolor="rgb(255,255,255"
    )
    axx$title <- "Time"
    axy$title <- "u"
    axz$title <- paste("Variable",var)
    Pl <- plot_ly(z = z0, x=time, y = u) %>%
      layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz))%>%
      add_surface(showscale=FALSE)
    print(Pl)
  } else stop("The type is not valid.")
}

