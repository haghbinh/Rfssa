#--------------------------------------------------------------
#' Functional Singular Spectrum Analysis Forecasting Bootstrap Prediction Interval
#'
#' This function calculates the bootstrap prediction interval for functional singular spectrum analysis (FSSA) forecasting predictions of univariate functional time series (\code{\link{fts}}) observed over a one-dimensional domain.
#' @return A list of numeric vectors where the first vector contains the discretely sampled point forecast, the second contains the discretely sampled lower bound of the prediction interval, and the third contains the discretely sampled upper bound of the prediction interval.
#' @param Y An object of class \code{\link{fts}}.
#' @param O A positive integer specifying the training set size.
#' @param L A positive integer specifying the window length.
#' @param ntriples The number of eigentriples to use to perform the forecasts.
#' @param Bt A positive integer specifying the number of bootstrap samples.
#' @param h An integer that specifies the forecast horizon.
#' @param alpha A double between zero and one specifying the significance level.
#' @param method A character string specifying the type of forecasting to perform either \code{"recurrent"} or \code{"vector"}.
#' @param tol A double specifying the amount of tolerated error in the approximation of the matrix that corresponds with the operator formed using a Neumann series leveraged in both forecasting algorithms.
#' @importFrom rainbow fts
#' @importFrom ftsa quantile.fts
#' @examples
#' \dontrun{
#' require(Rfssa)
#' if(!("ggplot2" %in% base::rownames(utils::installed.packages()))==TRUE){utils::install.packages("ggplot2")}
#' require(ggplot2)
#' load_github_data("https://github.com/haghbinh/Rfssa/blob/master/data/Callcenter.RData")
#' ## Define functional objects
#' D <- matrix(sqrt(Callcenter$calls), nrow = 240)
#' N <- ncol(D)
#' time <- substr(seq(ISOdate(1999, 1, 1), ISOdate(1999, 12, 31), by = "day"),1,10)
#' K <- nrow(D)
#' u <- seq(0, K, length.out = K)
#' d <- 22
#' ## Define functional time series
#' Y <- Rfssa::fts(list(D), list(list(d, "bspline")), list(u),time)
#' plot(Y, mains = c("Call Center Data Line Plot"),
#'      xlabels = "Time (6 minutes aggregated)",
#'      ylabels = "Sqrt of Call Numbers",type="line",
#'      xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),xticklocs =
#'        list(c(1,60,120,180,240)))
#' # Calculate prediction intervals
#'
#' pred_interval <- fpredinterval(Y = Y, O = 310, L = 28, ntriples = 7, Bt = 10000, h=3, alpha = 0.05, method = "recurrent")
#'
#' # Plot the forecast and prediction interval using ggplot
#'
#' fssa_forecast <- pred_interval$forecast
#' fssa_forecast_lower <- pred_interval$lower
#' fssa_forecast_upper <- pred_interval$upper
#' ggplot(data = data.frame(cbind(fssa_forecast,fssa_forecast_lower,fssa_forecast_upper)),aes(1:K,fssa_forecast))+theme_bw()+
#' geom_ribbon(aes(ymin=fssa_forecast_lower,ymax=fssa_forecast_upper),fill="grey90")+
#'   geom_line(aes(1:K,fssa_forecast_upper),color="black",size=0.5)+
#'   geom_line(aes(1:K,fssa_forecast_lower),color="black",size=0.5)+
#'   geom_line(aes(1:K,fssa_forecast),color="red2",size=1)+
#'   scale_x_continuous(name = "Time (6 minutes aggregated)", breaks=c(1,60,120,180,240),
#'                      labels = c("00:00","06:00","12:00","18:00","24:00"),)+
#'   scale_y_continuous(name="Sqrt of Call Numbers")+
#'   theme(axis.text.x=element_text(size=15),
#'         axis.title.x=element_text(size=18))+
#'   theme(axis.text.y=element_text(size=15),
#'         axis.title.y=element_text(size=18))+ggtitle("Prediction Intervals for Jan. 3, 2000")+
#'   theme(plot.title = element_text(hjust = 0.5,size=20))
#' }
#'
#' @export

fpredinterval <- function(Y, O, L, ntriples, Bt, h=1, alpha = 0.05, method = "recurrent", tol = 10^-3){
  M <- O+h
  g <- 1:ntriples
  basis <- Y@B[[1]]
  grid <- Y@grid[[1]]
  N <- ncol(Y@C[[1]])
  D <- basis%*%Y@C[[1]]
  E <- sapply(X = 1:(N-M), function(i){D[,(M+i)]-(basis%*%Rfssa::fforecast(U=Rfssa::fssa(Rfssa::fts(X = list(D[,i:(M+i-h)]),B = list(basis),grid = list(grid)),L=L),groups = list(g),h = h,method = method,tol = tol)[[1]]@C[[1]][,h])})
  E_B <- matrix(data = NA, nrow = length(grid), ncol = Bt)
  for(j in 1:Bt){E_B[,j] <- E[,sample(1:(N-M),1)]}
  colnames(E_B)=as.character(1:Bt)
  q_fts <- rainbow::fts(x = 1:length(grid), y = E_B)
  quants <- ftsa::quantile.fts(q_fts,probs = seq(0,1,0.01))
  upper_half <- 1-alpha/2
  lower_half <- alpha/2
  if((100*alpha)%%2==0){
    lower <- quants[,paste0(as.character(100*(lower_half)),"%")]
    upper <- quants[,paste0(as.character(100*(upper_half)),"%")]
  }else{
    lower <- (quants[,paste0(as.character(100*(lower_half+0.005)),"%")] + quants[,paste0(as.character(100*(lower_half-0.005)),"%")])/2
    upper <- (quants[,paste0(as.character(100*(upper_half+0.005)),"%")] + quants[,paste0(as.character(100*(upper_half-0.005)),"%")])/2
  }
  fssa_forecast <- basis%*%Rfssa::fforecast(Rfssa::fssa(Y,L=L),groups = list(g),h = h,method = method)[[1]]@C[[1]][,h]
  fssa_forecast_lower <- fssa_forecast+lower
  fssa_forecast_upper <- fssa_forecast+upper

  return(list(forecast=fssa_forecast,lower=fssa_forecast_lower,upper=fssa_forecast_upper))

}
