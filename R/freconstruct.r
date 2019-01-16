#--------------------------------------------------------------
#' Reconstruction sage of FSSA
#'
#' This is a function for the reconstructing functional time series from
#'functional singular spectrum objects (including Grouping and
#' Hankelization steps). The output is a list of functional time series corresponds to each group.
#' @return a named list of reconstructed functional time series in each groups and
#' a numeric vector of eigenvalues.
#' @param U an object of class \code{\link{fssa}}
#' @param group a list of numeric vectors, each vector includes indices of such elementary components
#' of a group used for reconstruction.
#' @seealso \code{\link{fssa}}

#' @importFrom fda fd

#' @export
freconstruct <- function(U, group = as.list(1L:10L)) {
  N <- U$N
  Y <- U$Y
  d <- nrow(U[[1]]$coefs)
  L <- U$L
  K <- N - L + 1L
  basis <- U[[1]]$basis
  m <- length(group)
  basis <- U[[1]]$basis
  out <- list()
  for (i in 1L:m) {
    Cx <- matrix(NA, nrow = d, ncol = N)
    g <- group[[i]]
    S <- 0L
    for (j in 1L:length(g)) S <- S +
      fproj(U, g[j], d, K, L, Y)
    S <- fH(S, d)
    Cx[, 1L:L] <- S[, 1L, ]
    Cx[, L:N] <- S[, ,L]
    out[[i]] <- fd(Cx, basis)
  }
  out$values <- sqrt(U$values)
  return(out)
}
