#--------------------------------------------------------------
#' Reconstruction Stage of Functional Singular Spectrum Analysis
#'
#' This is a function for reconstructing functional time series (\code{\link{fts}}) objects from functional singular spectrum analysis
#' (\code{\link{fssa}}) objects (including Grouping and
#' Hankelization steps). The function performs the reconstruction step for univariate functional singular spectrum analysis (ufssa) or multivariate functional singular spectrum analysis (mfssa)
#' depending on whether or not the input is an \code{\link{fssa}} object from ufssa or mfssa.
#' @return a named list of reconstructed functional time series in each groups and
#' a numeric vector of eigenvalues
#' @param U an object of class \code{\link{fssa}}
#' @param group a list of numeric vectors, each vector includes indices of elementary components
#' of a group used for reconstruction
#' @seealso \code{\link{fssa}} for an example on how to run this function starting from \code{\link{fssa}} objects

#' @export
freconstruct <- function(U, group = as.list(1L:10L)) {
  if(is.fd(U[[1]])) out <- ufreconstruct(U,group) else out <- mfreconstruct(U,group)
  return(out)
}
