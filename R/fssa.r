#' Functional Singular Spectrum Analysis
#'
#' This is a function which performs the decomposition (including embedding
#'  and  functional SVD steps) stage for univariate functional singular spectrum analysis (ufssa)
#'  or multivariate functional singular spectrum analysis (mfssa). The algorithm (ufssa or mfssa) is chosen based on
#'  whether the supplied input is a univariate or
#'  multivariate functional time series (\code{\link{fts}}) object. The \code{type} parameter can also be set to \code{"mfssa"} if the user wishes to perform ufssa of a univariate \code{\link{fts}} object using the mfssa code. Also note that the variables of the \code{\link{fts}} maybe observed over different dimensional domains where the maximum dimension currently supported is two.
#' @return An object of class \code{fssa}, which is a list of functional objects and the following components:
#' \item{values}{A numeric vector of eigenvalues.}
#' \item{L}{The specified window length.}
#' \item{N}{The length of the functional time series.}
#' \item{Y}{The original functional time series.}
#' @param Y An object of class \code{\link{fts}}.
#' @param L A positive integer giving the window length.
#' @param ntriples A positive integer specifying the number of eigentriples to calculate in the decomposition.
#' @param type A string indicating which type of fssa to perform. Use \code{type="ufssa"} to perform univariate fssa (default for univariate fts). Use \code{type="mfssa"} to perform multivariate fssa (default for multivariate fts).
#' @importFrom RSpectra eigs
#' @examples
#' \dontrun{
#' ## Univariate FSSA Example on Callcenter data
#' require(Rfssa)
#' load_github_data("https://github.com/haghbinh/Rfssa/blob/master/data/Callcenter.RData")
#' ## Define functional objects
#' D <- matrix(sqrt(Callcenter$calls), nrow = 240)
#' N <- ncol(D)
#' time <- substr(seq(ISOdate(1999, 1, 1), ISOdate(1999, 12, 31), by = "day"),1,10)
#' K <- nrow(D)
#' u <- seq(0, K, length.out = K)
#' d <- 22 # Optimal Number of basis elements
#' ## Define functional time series
#' Y <- Rfssa::fts(list(D), list(list(d, "bspline")), list(u),time)
#' Y
#' plot(Y, mains = c("Call Center Data Line Plot"),
#' xlabels = "Time (6 minutes aggregated)",
#' ylabels = "Sqrt of Call Numbers",type="line",
#' xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),xticklocs =
#'  list(c(1,60,120,180,240)))
#' ## Univariate functional singular spectrum analysis
#' L <- 28
#' U <- fssa(Y, L)
#' plot(U, d = 13)
#' plot(U, d = 9, type = "lheats")
#' plot(U, d = 9, type = "lcurves")
#' plot(U, d = 9, type = "vectors")
#' plot(U, d = 10, type = "periodogram")
#' plot(U, d = 10, type = "paired")
#' plot(U, d = 10, type = "wcor")
#' gr <- list(1, 2:3, 4:5, 6:7, 1:7)
#' Q <- freconstruct(U, gr)
#' plot(Q[[1]], mains = "Call Center Mean Component",
#' xlabels = "Time (6 minutes aggregated)",
#'      ylabels = "Sqrt of Call Numbers",type="line",
#'      xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),xticklocs =
#'        list(c(1,60,120,180,240)))
#' plot(Q[[2]], mains = "Call Center First Periodic Component",
#' xlabels = "Time (6 minutes aggregated)",
#'      ylabels = "Sqrt of Call Numbers",type="line",
#'      xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),xticklocs =
#'        list(c(1,60,120,180,240)))
#' plot(Q[[3]], mains = "Call Center Second Periodic Component",
#' xlabels = "Time (6 minutes aggregated)",
#'      ylabels = "Sqrt of Call Numbers",type="line",
#'      xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),xticklocs =
#'        list(c(1,60,120,180,240)))
#' plot(Q[[4]], mains = "Call Center Third Periodic Component",
#' xlabels = "Time (6 minutes aggregated)",
#'      ylabels = "Sqrt of Call Numbers",type="line",
#'      xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),xticklocs =
#'        list(c(1,60,120,180,240)))
#' plot(Y, mains = c("Call Center Data Line Plot"),
#' xlabels = "Time (6 minutes aggregated)",
#'      ylabels = "Sqrt of Call Numbers",type="line",
#'      xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),xticklocs =
#'        list(c(1,60,120,180,240)))
#' plot(Q[[5]], mains = "Call Center Extracted Signal",
#' xlabels = "Time (6 minutes aggregated)",
#'      ylabels = "Sqrt of Call Numbers",type="line",
#'      xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),xticklocs =
#'        list(c(1,60,120,180,240)))
#'
#' ## Other visualization types for object of class "fts":
#'
#' plot(Q[[1]],xlabels = "Time (6 minutes aggregated)",
#'      zlabels = "Sqrt of Call Numbers",type="3Dsurface", tlabels = "Date",
#'      xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),xticklocs =
#'        list(c(1,60,120,180,240)))
#' plot(Q[[2]], mains = "Call Center First Periodic Component",
#' xlabels = "Time (6 minutes aggregated)",
#'      zlabels = "Sqrt of Call Numbers",type="heatmap",
#'      tlabels = "Date",
#'      xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),xticklocs =
#'        list(c(1,60,120,180,240)))
#' plot(Q[[3]],xlabels = "Time (6 minutes aggregated)",
#'      zlabels = "Sqrt of Call Numbers",type="3Dline",
#'      tlabels = "Date",
#'      xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),xticklocs =
#'        list(c(1,60,120,180,240)))
#'
#' ## Multivariate FSSA Example on bivariate intraday
#' ## temperature curves and smoothed images of vegetation
#' require(Rfssa)
#' load_github_data("https://github.com/haghbinh/Rfssa/blob/master/data/Montana.RData")
#' Temp <- Montana$Temp
#' NDVI <- Montana$NDVI
#' d_temp <- 11
#' d_NDVI <- 13
#' ## Define functional time series
#' Y <- Rfssa::fts(
#'   list(Temp / sd(Temp), NDVI), list(
#'     list(d_temp, "bspline"),
#'     list(d_NDVI, d_NDVI, "bspline", "bspline")
#'   ),
#'   list(c(0, 23), list(c(1, 33), c(1, 33))),
#' colnames(Temp))
#' # Plot the first 100 observations
#' plot(Y[1:100],
#'      xlabels = c("Time", "Longitude"),
#'      ylabels = c("Normalized Temperature (\u00B0C)", "Latitude"),
#'      zlabels = c("", "NDVI"),
#'      mains = c("Temperature Curves", "NDVI Images"), color_palette = "RdYlGn",
#'      xticklabels = list(c("00:00","06:00","12:00","18:00","24:00"),
#'      c("113.40\u00B0 W", "113.30\u00B0 W")),xticklocs =
#'        list(c(1,6,12,18,24),c(1,33)),
#'        yticklabels = list(NA,c("48.70\u00B0 N", "48.77\u00B0 N")),yticklocs =
#'        list(NA,c(1,33))
#' )
#' plot(Y, types = c("3Dline", "heatmap"), vars = c(1, 1))
#' plot(Y, types = "heatmap", vars = 2)
#' plot(Y, vars = c(2, 1))
#' L <- 45
#' ## Multivariate functional singular spectrum analysis
#' U <- fssa(Y, L)
#' plot(U, type = "values", d = 10)
#' plot(U, type = "vectors", d = 4)
#' plot(U, type = "lheats", d = 4)
#' plot(U, type = "lcurves", d = 4, vars = c(1))
#' plot(U, type = "paired", d = 6)
#' plot(U, type = "wcor", d = 10)
#' plot(U, type = "periodogram", d = 4)
#' # Reconstruction of multivariate fts observed over different dimensional domains
#' Q <- freconstruct(U = U, groups = list(c(1), c(2:3), c(4)))
#' # Plotting reconstructions to show accuracy
#' plot(Q[[1]],
#'      xlabels = c("Time", "Longitude"),
#'      ylabels = c("Normalized Temperature (\u00B0C)", "Latitude"),
#'      zlabels = c("", "NDVI"),
#'      mains = c("Temperature Curves Mean", "NDVI Images Mean"), color_palette = "RdYlGn",
#'      xticklabels = list(c("00:00","06:00","12:00","18:00","24:00"),
#'      c("113.40\u00B0 W", "113.30\u00B0 W")),xticklocs =
#'        list(c(1,6,12,18,24),c(1,33)),
#'        yticklabels = list(NA,c("48.70\u00B0 N", "48.77\u00B0 N")),yticklocs =
#'        list(NA,c(1,33))) # mean
#' plot(Q[[2]],
#'      xlabels = c("Time", "Longitude"),
#'      ylabels = c("Normalized Temperature (\u00B0C)", "Latitude"),
#'      zlabels = c("", "NDVI"),
#'      mains = c("Temperature Curves Periodic", "NDVI Images Periodic"), color_palette = "RdYlGn",
#'      xticklabels = list(c("00:00","06:00","12:00","18:00","24:00"),
#'      c("113.40\u00B0 W", "113.30\u00B0 W")),xticklocs =
#'        list(c(1,6,12,18,24),c(1,33)),
#'        yticklabels = list(NA,c("48.70\u00B0 N", "48.77\u00B0 N")),yticklocs =
#'        list(NA,c(1,33))) # periodic
#' plot(Q[[3]],
#'      xlabels = c("Time", "Longitude"),
#'      ylabels = c("Normalized Temperature (\u00B0C)", "Latitude"),
#'      zlabels = c("", "NDVI"),
#'      mains = c("Temperature Curves Trend", "NDVI Images Trend"), color_palette = "RdYlGn",
#'      xticklabels = list(c("00:00","06:00","12:00","18:00","24:00"),
#'      c("113.40\u00B0 W", "113.30\u00B0 W")),xticklocs =
#'        list(c(1,6,12,18,24),c(1,33)),
#'        yticklabels = list(NA,c("48.70\u00B0 N", "48.77\u00B0 N")),yticklocs =
#'        list(NA,c(1,33))) # trend
#' }
#' @useDynLib Rfssa
#' @export
fssa <- function(Y, L = NA, ntriples = 20, type = "ufssa") {
  N <- Y$N
  if (is.na(L)) L <- floor(N / 2L)
  if (ntriples > L){
    ntriples = L
    warning("\"ntriples\" must be less than or equal to \"L\". Setting \"ntriples\" = \"L\"")
  }
  p <- length(Y$dimSupp)
  if (p == 1 && type == "ufssa") {
    out <- ufssa(Y, L, ntriples)
  } else if (p > 1 || type == "mfssa") {
    out <- mfssa(Y, L, ntriples)
  } else {
    stop("Error in type or dimension.")
  }
  class(out) <- "fssa"
  return(out)
}
