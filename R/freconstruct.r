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
  if(is.fd(U[[1]])) out <- ufreconstruct(U,group) else out <- mfreconstruct(U,group)
  return(out)
}
