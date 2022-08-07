#--------------------------------------------------------------
#' Reconstruction Stage of Functional Singular Spectrum Analysis
#'
#' This is a function for reconstructing univariate or multivariate functional time series (\code{\link{fts}}) objects from functional singular spectrum analysis
#' (\code{\link{fssa}}) objects (including Grouping and
#' Hankelization steps). The function performs the reconstruction step for univariate functional singular spectrum analysis (ufssa) or multivariate functional singular spectrum analysis (mfssa)
#' depending on whether or not the input is an \code{\link{fssa}} object from ufssa or mfssa.
#' @return A named list of objects of class \code{\link{fts}} that are reconstructed as according to the specified groups and
#' a numeric vector of eigenvalues.
#' @param U An object of class \code{\link{fssa}}.
#' @param groups A list of numeric vectors, each vector includes indices of elementary components.
#' of a group used for reconstruction.
#' @note Refer to \code{\link{fssa}} for an example on how to run this function starting from \code{\link{fssa}} objects.
#' @seealso \code{\link{fssa}}, \code{\link{fts}},
#' @export
freconstruct <- function(U, groups = as.list(1L:10L)) {
  if (class(U[[1]])[[1]] != "list") out <- ufreconstruct(U, groups) else out <- mfreconstruct(U, groups)
  return(out)
}
