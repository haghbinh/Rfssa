#--------------------------------------------------------------
#' FSSA Forecasting Bootstrap Prediction Interval
#'
#' Calculate the bootstrap prediction interval for functional singular
#' spectrum analysis (FSSA) forecasting predictions of univariate functional
#' time series (\code{\link{funts}}) observed over a one-dimensional domain.
#' @param Y an object of class \code{\link{funts}}.
#' @param O a positive integer specifying the training set size.
#' @param L a positive integer specifying the window length.
#' @param ntriples the number of eigentriples to use for forecasts.
#' @param Bt a positive integer specifying the number of bootstrap samples.
#' @param h an integer specifying the forecast horizon.
#' @param alpha a double (0 < alpha < 1) specifying the significance level.
#' @param method a character string: "recurrent" or "vector" forecasting.
#' @param tol a double specifying tolerated error in the approximation.
#' @return a list of numeric vectors: point forecast, lower, and upper bounds.


#' @importFrom rainbow fts
#' @importFrom ftsa quantile.fts
#' @examples
#' \dontrun{
#' data("Callcenter")
#' pred_interval <- fpredinterval(
#'   Y = Callcenter, O = 310,
#'   L = 28, ntriples = 7, Bt = 10000, h = 3
#' )
#'
#' # Plot the forecast and prediction interval using ggplot
#' df <- data.frame(
#'   x = 1:240,
#'   y = pred_interval$forecast,
#'   lower = pred_interval$lower,
#'   upper = pred_interval$upper
#' )
#' require(ggplot2)
#' # Create the ggplot
#' ggplot(df, aes(x = x, y = y)) +
#'   geom_line(size = 1.2) +
#'   scale_x_continuous(
#'     name = "Time",
#'     breaks = c(1, 60, 120, 180, 240),
#'     labels = c("00:00", "06:00", "12:00", "18:00", "24:00"),
#'   ) +
#'   scale_y_continuous(name = "Sqrt of Call Numbers") +
#'   ggtitle("Prediction Intervals for Jan. 3, 2000") +
#'   geom_ribbon(aes(ymin = lower, ymax = upper), fill = "darkolivegreen3", alpha = 0.3) +
#'   theme_minimal()
#' }
#'
#' @export

fpredinterval <- function(Y, O, L, ntriples, Bt, h = 1, alpha = 0.05, method = "recurrent", tol = 10^-3) {
  cat("Running, please wait...\n")
  N <- Y$N
  start_t <- Y$time[1]
  end_t <- Y$time[N]
  basisobj <- Y$basis
  argval <- Y$argval
  p <- length(Y$dimSupp)
  if (p == 1) {
    basisobj <- basisobj[[1]]
    argval <- argval[[1]]
  }
  M <- O + h - 1
  g <- 1:ntriples
  basis <- Y$B_mat[[1]]
  grid <- Y$argval[[1]]
  N <- ncol(Y$coefs[[1]])
  D <- basis %*% Y$coefs[[1]]
  E <- sapply(X = 1:(N - M), function(i) {
    x_funts <- funts(X = D[, i:(M + i - h)], basisobj = basisobj, argval = argval, start = start_t, end = end_t)
    if (p == 1) {
      U <- ufssa(x_funts, L = L, 20)
      fore <- ufforecast(U, groups = list(g), h = h, method = method, tol = tol)
    } else {
      U <- mfssa(x_funts, L = L, 20)
      fore <- mfforecast(U, groups = list(g), h = h, method = method, tol = tol)
    }
    D[, (M + i)] - (basis %*% fore[[1]]$coefs[[1]][, h])
  })
  E_B <- matrix(data = NA, nrow = length(grid), ncol = Bt)
  for (j in 1:Bt) {
    E_B[, j] <- E[, sample(1:(N - M), 1)]
  }
  colnames(E_B) <- as.character(1:Bt)
  q_fts <- rainbow::fts(x = 1:length(grid), y = E_B)
  quants <- ftsa::quantile.fts(q_fts, probs = seq(0, 1, 0.01))
  upper_half <- 1 - alpha / 2
  lower_half <- alpha / 2
  if ((100 * alpha) %% 2 == 0) {
    lower <- quants[, paste0(as.character(100 * (lower_half)), "%")]
    upper <- quants[, paste0(as.character(100 * (upper_half)), "%")]
  } else {
    lower <- (quants[, paste0(as.character(100 * (lower_half + 0.005)), "%")] + quants[, paste0(as.character(100 * (lower_half - 0.005)), "%")]) / 2
    upper <- (quants[, paste0(as.character(100 * (upper_half + 0.005)), "%")] + quants[, paste0(as.character(100 * (upper_half - 0.005)), "%")]) / 2
  }
  if (p == 1) {
    U <- ufssa(Y, L = L, 20)
    fore <- ufforecast(U, groups = list(g), h = h, method = method, tol = tol)
  } else {
    U <- mfssa(Y, L = L, 20)
    fore <- mfforecast(U, groups = list(g), h = h, method = method, tol = tol)
  }
  fssa_forecast <- basis %*% fore[[1]]$coefs[[1]][, h]
  fssa_forecast_lower <- fssa_forecast + lower
  fssa_forecast_upper <- fssa_forecast + upper
  cat("Done.\n")
  return(list(forecast = fssa_forecast, lower = fssa_forecast_lower, upper = fssa_forecast_upper))
}
