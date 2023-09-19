#' Weighted-Correlations Plot
#'
#' This function generates a plot displaying the weighted-correlation (w-correlation) matrix of functional time series (\code{\link{fts}}) objects reconstructed from functional singular spectrum analysis (\code{\link{fssa}}) objects.
#'
#' @param W a w-correlation matrix.
#' @param cuts an integer specifying the number of levels for dividing the range of w-correlation values.
#' @param main a character string for the main title of the plot.
#'
#' @note For an example of how to use this function with a w-correlation matrix, refer to \code{\link{fwcor}}.
#'
#' @seealso \code{\link{fssa}}, \code{\link{freconstruct}}, \code{\link{funts}}, \code{\link{fwcor}}
#' @export
wplot <- function(W, cuts = 20, main = NA) {
  at <- pretty(c(0, 1), n = cuts)
  d <- nrow(W)
  W0 <- abs(W)
  a <- min(W0)
  b <- max(W0 - diag(1, d))
  s <- stats::sd(W0 - diag(1, d))
  diag(W0) <- min(1, b + 3 * s)
  xylabels <- paste0("F", 1:d)
  if (is.na(main)) main <- "W-correlation matrix"
  p1 <- lattice::levelplot(1 - W0,
    xlab = "", at = at,
    ylab = "", colorkey = NULL,
    main = list(main, cex = 2),
    scales = list(
      x = list(
        at = 1:d,
        lab = xylabels, cex = 0.9
      ),
      y = list(
        at = 1:d,
        lab = xylabels, cex = 0.9
      )
    ),
    col.regions = grDevices::gray(seq(0, 1, length = 100))
  )
  graphics::plot(p1)
}
