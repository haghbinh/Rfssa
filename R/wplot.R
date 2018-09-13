#' W-correlation matrix plot
#'
#' W-correlation matrix of the single components.
#' @param W a matrix of W-correlation values.
#' @export
wplot <- function(W) {
  d <- nrow(W)
  W0 <- abs(W)
  a <- min(W0)
  b <- max(W0 - diag(1, d))
  s <- stats::sd(W0 - diag(1, d))
  diag(W0) <- min(1, b + 3 * s)
  xylabels <- paste0("F", 1:d)
  p1 <- lattice::levelplot(1 - W0, xlab = "",
                           ylab = "", colorkey = NULL,
                           main = paste("W-correlation matrix"),
                           scales = list(x = list(at = 1:d,
                                                  lab = xylabels),
                                         y = list(at = 1:d,
                                                  lab = xylabels)),
                           col.regions = grDevices::gray(seq(0,                                                                                                    1, length = 100)))
  graphics::plot(p1)
}
