#--------------------------------------------------------------
#' Functional Singular Spectrum Analysis Recurrent Forecasting and Vector Forecasting
#'
#' This function performs functional singular spectrum analysis (FSSA) recurrent forecasting (FSSA R-forecasting) or vector forecasting (FSSA V-forecasting) of univariate or multivariate functional time series (\code{\link{fts}}) observed over a one-dimensional domain.
#' @return A list of objects of class \code{\link{fts}} where each fts corresponds to a forecasted group.
#' @param U An object of class \code{\link{fssa}} that holds the decomposition.
#' @param groups A list of numeric vectors where each vector includes indices of elementary components of a group used for reconstruction and forecasting.
#' @param h An integer that specifies the forecast horizon.
#' @param method A character string specifying the type of forecasting to perform either \code{"recurrent"} or \code{"vector"}.
#' @param tol A double specifying the amount of tolerated error in the approximation of the matrix that corresponds with the operator formed using a Neumann series leveraged in both forecasting algorithms; see Trinka et. al. (2021) for more details.
#' @examples
#' \dontrun{
#' # FSSA Forecasting
#' require(Rfssa)
#' load_github_data("https://github.com/haghbinh/Rfssa/blob/master/data/Callcenter.RData")
#' ## Define functional objects
#' D <- matrix(sqrt(Callcenter$calls), nrow = 240)
#' N <- ncol(D)
#' time <- seq(ISOdate(1999, 1, 1), ISOdate(1999, 12, 31), by = "day")
#' K <- nrow(D)
#' u <- seq(0, K, length.out = K)
#' d <- 22
#' ## Define functional time series
#' Y <- fts(list(D), list(list(d, "bspline")), list(u))
#' plot(Y, mains = "Call Center Data")
#' ## Perfom FSSA decomposition
#' L <- 28
#' U <- fssa(Y, L)
#' groups <- list(1:7, 1, 2:3, 4:5, 6:7)
#' ## Perform FSSA R-forecast and FSSA V-forecast
#' pr_R <- fforecast(U = U, groups = groups, h = 30, method = "recurrent", tol = 10^-3)
#' plot(pr_R[[1]], mains = "Recurrent Forecast Group 1")
#' plot(pr_R[[2]], mains = "Recurrent Forecast Group 2")
#' plot(pr_R[[3]], mains = "Recurrent Forecast Group 3")
#' plot(pr_R[[4]], mains = "Recurrent Forecast Group 4")
#' plot(pr_R[[5]], mains = "Recurrent Forecast Group 5")
#'
#' pr_V <- fforecast(U = U, groups = groups, h = 30, method = "vector", tol = 10^-3)
#' plot(pr_V[[1]], mains = "Vector Forecast Group 1")
#' plot(pr_V[[2]], mains = "Vector Forecast Group 2")
#' plot(pr_V[[3]], mains = "Vector Forecast Group 3")
#' plot(pr_V[[4]], mains = "Vector Forecast Group 4")
#' plot(pr_V[[5]], mains = "Vector Forecast Group 5")
#'
#' # MFSSA Forecasting
#' require(Rfssa)
#' load_github_data("https://github.com/haghbinh/Rfssa/blob/master/data/Jambi.RData")
#' ## Raw image data
#' NDVI <- Jambi$NDVI
#' EVI <- Jambi$EVI
#' time <- Jambi$Date
#' ## Kernel density estimation of pixel intensity
#' D0_NDVI <- matrix(NA, nrow = 512, ncol = 448)
#' D0_EVI <- matrix(NA, nrow = 512, ncol = 448)
#' for (i in 1:448) {
#'   D0_NDVI[, i] <- density(NDVI[, , i], from = 0, to = 1)$y
#'   D0_EVI[, i] <- density(EVI[, , i], from = 0, to = 1)$y
#' }
#' ## Define functional objects
#' d <- 11
#' Y <- fts(list(D0_NDVI, D0_EVI), list(list(d, "bspline"), list(d + 4, "fourier")), list(c(0, 1), c(0, 1)))
#' plot(Y)
#' U <- fssa(Y = Y, L = 45)
#' groups <- list(c(1:4), c(1), c(2:3), c(4))
#' pr_R <- fforecast(U = U, groups = groups, h = 1, method = "recurrent")
#' plot(pr_R[[1]])
#' plot(pr_R[[2]])
#' plot(pr_R[[3]])
#' plot(pr_R[[4]])
#'
#' pr_V <- fforecast(U = U, groups = groups, h = 1, method = "vector")
#' plot(pr_V[[1]])
#' plot(pr_V[[2]])
#' plot(pr_V[[3]])
#' plot(pr_V[[4]])
#' }
#'
#' @export
fforecast <- function(U, groups = list(c(1)), h = 1, method = "recurrent", tol = 10^-3) {
  for (j in 1:length(U$Y@C)) {
    if (ncol(U$Y@grid[[j]]) > 1) {
      stop("Current forecasting routines only support fts whose variables are observed over one-dimensional domains. Forecasting of fts variables whose domains have dimension greater than one is under development.")
    }
  }

  if (length(U$Y@C) == 1) {
    out <- ufforecast(U = U, groups = groups, h = h, method = method, tol = tol)
  } else {
    out <- mfforecast(U = U, groups = groups, h = h, method = method, tol = tol)
  }

  return(out)
}
