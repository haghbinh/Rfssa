#' Correlation for Functional Time Series Objects
#'
#' This function finds the correlation between univariate or multivariate functional time series (\code{\link{fts}}) objects.
#' @return A scalar that is the correlation between \code{\link{fts}} objects.
#' @param Y1 An object of class \code{\link{fts}}.
#' @param Y2 An object of class \code{\link{fts}}.
#'
#'
#' @seealso \code{\link{fts}}
#'
#' @examples
#' \dontrun{
#' require(Rfssa)
#' load_github_data("https://github.com/haghbinh/Rfssa/blob/master/data/Jambi.RData")
#' ## Raw image data
#' NDVI <- Jambi$NDVI
#' EVI <- Jambi$EVI
#' ## Kernel density estimation of pixel intensity
#' D0_NDVI <- matrix(NA, nrow = 512, ncol = 448)
#' D0_EVI <- matrix(NA, nrow = 512, ncol = 448)
#' for (i in 1:448) {
#'   D0_NDVI[, i] <- density(NDVI[, , i], from = 0, to = 1)$y
#'   D0_EVI[, i] <- density(EVI[, , i], from = 0, to = 1)$y
#' }
#' d <- 11
#' u <- seq(0, 1, length.out = 512)
#' Y_1 <- Rfssa::fts(list(D0_NDVI), list(list(d, "bspline")), list(u))
#' Y_2 <- Rfssa::fts(list(D0_EVI), list(list(d, "bspline")), list(u))
#' out <- cor.fts(Y_1, Y_2)
#' print(out)
#' }
#'
#' @export
cor.fts <- function(Y1, Y2) {
  if (ncol(Y1@C[[1]]) != ncol(Y2@C[[1]])) {
    stop("Functional time series are of different length")
  }
  if (length(Y1) != length(Y2)) {
    stop("Functional time series have different number of covariates")
  }
  N <- ncol(Y1@C[[1]])
  w <- matrix(nrow = 1, ncol = N, data = 1)
  p <- length(Y1)
  G <- list()
  Y1_list <- list()
  Y2_list <- list()
  for (i in 1:p) {
    if (ncol(Y1@grid[[i]]) == 1) {
      G[[i]] <- t(onedG(A = Y1@B[[i]], B = Y1@B[[i]], grid = Y1@grid[[i]]))
    } else {
      G[[i]] <- t(twodG(A = Y1@B[[i]], B = Y1@B[[i]], grid = Y1@grid[[i]]))
    }
    Y1_list[[i]] <- Y1@C[[i]]
    Y2_list[[i]] <- Y2@C[[i]]
  }

  wcor <- mwinprod(Y1_list, Y2_list, w, G, p) / sqrt(mwinprod(Y1_list, Y1_list, w, G, p) * mwinprod(Y2_list, Y2_list, w, G, p))

  return(wcor)
}
