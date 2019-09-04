#' Rfssa: A package for functional singular spectrum analysis.
#'
#' The Rfssa package provides the collections of necessary
#' functions to implement Functional Singular Spectrum Analysis (FSSA)
#' for analysing Functional Time Series (FTS). FSSA is a novel
#' non-parametric method to perform decomposition and reconstruction of FTS.
#'
#'
#'@details
#'Typically the use of the package starts with the decomposition
#'of the functional time series using \code{\link{fssa}}.
#'Then,  a suitable grouping of the elementary time series is required.
#' This can be done heuristically, for example, via looking at
#' the plots of the decomposition (\code{\link[=plot.fssa]{plot}}).
#'  Alternatively, one can examine the so-called w-correlation matrix
#'  (\code{\link{fwcor}}).
#'  Next step includes the reconstruction of the time-series using
#'  the selected grouping (\code{\link{freconstruct}}).
#'
#'
#'@seealso
#'  \code{\link{fssa}}, \code{\link{freconstruct}},
#'  \code{\link{fwcor}}, \code{\link{wplot}}, \code{\link{plot.fts}}
#'
#'
#'@references
#'   Haghbin H., Najibi, S.M., Mahmoudvand R., Trinka J., Maadooliat M. (2019).
#'   Functional singular spectrum Analysis. Manuscript submitted for publication.

#' @docType package
#' @name Rfssa
NULL
