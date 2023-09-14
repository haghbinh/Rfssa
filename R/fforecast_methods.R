# =======================================================================
#'
#' Plot method for FSSA Forecast (fforecast) class
#'
#' Create visualizations of FSSA Forecast (fforecast) class. This function supports
#' plotting `fforecast` data with one-dimensional or two-dimensional domains.
#'
#' @param fobj an object of class \code{fforecast}.
#' @param group an integer specify the group index for plot.
#' @param npts number of grid points for the plots.
#' @param obs observation number (for two-dimensional domains).
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param main main title for the plot.
#' @param type type of plot ("l" for line, "p" for points, etc.).
#' @param lty line type (1 for solid, 2 for dashed, etc.).
#' @param ... additional graphical parameters passed to plotting functions.
#'
#' @seealso \code{\link{fforecast}}
#'
#' @export
plot.fforecast <- function(fobj, group = 1, npts = 100, obs = 1, xlab = NULL, ylab = NULL, main = NULL, type = "l", lty = 1, ...) {
  obj <- fobj[[group]]
  plot(obj, npts, obs, xlab, ylab, main, type, lty, ...)
}
