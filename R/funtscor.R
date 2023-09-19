#' Correlation for Functional Time Series (funts) Objects
#'
#' Find the correlation between univariate or multivariate functional time series
#' (\code{\link{funts}}) objects.
#'
#' @param Y1 an object of class \code{\link{funts}}.
#' @param Y2 an object of class \code{\link{funts}}.
#'
#' @return A scalar representing the correlation between \code{\link{funts}} objects.
#'
#' @seealso \code{\link{funts}}
#'
#' @examples
#' \dontrun{
#' require(Rfssa)
#' data("Montana")
#' ## Raw image data
#' NDVI <- Montana[, 1]
#' EVI <- Montana[, 2]
#' cor.funts(NDVI, EVI)
#' }
#'
#' @export
cor.funts <- function(Y1, Y2) {
  if (ncol(Y1$coefs[[1]]) != ncol(Y2$coefs[[1]])) {
    stop("Functional time series have different lengths")
  }
  if (length(Y1) != length(Y2)) {
    stop("Functional time series have different numbers of covariates")
  }
  N <- ncol(Y1$coefs[[1]])
  w <- matrix(nrow = 1, ncol = N, data = 1)
  p <- length(Y1$dimSupp)
  G <- list()
  Y1_list <- list()
  Y2_list <- list()

  for (i in 1:p) {
    if (Y1$dimSupp[[i]] == 1) {
      G[[i]] <- t(onedG(A = Y1$B_mat[[i]], B = Y1$B_mat[[i]], grid = Y1$argval[[i]]))
    } else {
      G[[i]] <- t(twodG(A = Y1$B_mat[[i]], B = Y1$B_mat[[i]], grid = Y1$argval[[i]]))
    }
    Y1_list[[i]] <- Y1$coefs[[i]]
    Y2_list[[i]] <- Y2$coefs[[i]]
  }

  wcor <- mwinprod(Y1_list, Y2_list, w, G, p) / sqrt(mwinprod(Y1_list, Y1_list, w, G, p) * mwinprod(Y2_list, Y2_list, w, G, p))

  return(wcor)
}
