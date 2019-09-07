#' Rfssa: A Package for Functional Singular Spectrum Analysis.
#'
#' The Rfssa package provides the collection of necessary
#' functions to implement functional singular spectrum analysis (FSSA)
#' for analysing functional time series (FTS) data. FSSA is a novel
#' non-parametric method to perform decomposition and reconstruction of FTS.
#'
#'
#'@details
#' The use of the package starts with the decomposition
#' of the functional time series using \code{\link{fssa}}.
#' Then a suitable grouping of the principal compomnents is required for reconstruction.
#' This can be done heuristically by looking at
#' the plots of the decomposition (\code{\link[=plot.fssa]{plot}}).
#' Alternatively, one can examine the  weighted correlation (w-correlation) matrix
#' (\code{\link{fwcor}}). The final step is the reconstruction of the
#' principal components into additive functional time series
#' whose sum approximates the original univariate or multivariate functional time series (\code{\link{freconstruct}}).
#'
#' This version of the Rfssa package includes updates to existing functions including \code{\link{fssa}}, \code{\link[=plot.fssa]{plot}}, \code{\link{wplot}}, and
#' \code{\link{freconstruct}}. Multivariate functional singular spectrum analysis (mfssa) was added to the package in \code{\link{fssa}} to allow
#' the user to perform embedding and decomposition of a multivariate FTS. The reconstruction stage in \code{\link{freconstruct}} was also updated
#'  to allow for reconstruction (including Hankelization) of multivariate FTS objects
#'  using multivariate FSSA objects that come from mfssa. Plotting options for FSSA objects in \code{\link[=plot.fssa]{plot}} were also updated
#'  so that the user can now plot left singular functions, right singular vectors, left singular function heat diagrams, and periodograms. FSSA plotting
#'  options also allow the user to specify which particular components they want to plot. For example, a user can specify that they want to see a paired-plot of only the
#'  third and fourth component. The 'meanvectors' and 'meanpaired' options were removed as these are satisfied with 'paired' and 'vectors' options. The 'efunctions'
#'  and 'efunctions2' options were also removed in lieu of the addition of the left singular function heat map option. The user can also specify the 'cuts' parameter in
#'  \code{\link{wplot}} to make visualization of the w-correlation matrix easier.
#'
#' This version of the Rfssa package also includes new functions
#' for converting functional data (FD) objects to FTS objects, arithmetic, indexing, correlation, and plotting of FTS data.
#' The user is able to convert an FD object to an FTS object using \code{\link{fts}}. The user can also perform addition,
#' subtraction, and multiplication of FTS objects with other FTS objects or FTS objects with scalars
#' by using '+', '-', and '*' respectively. The package also
#' allows for indexing of FTS objects by using '['. The user can also measure the unweighted correlation
#' between FTS objects by using \code{\link{cor.fts}}. The plotting of FTS objects can be performed using \code{\link[=plot.fts]{plot}}
#' which uses the plotly package for visualization .
#'
#' The package update also includes a new shiny app (\code{\link{launchApp}}) that can be used for demonstrations of univariate or multivariate FSSA
#' depending on the type that is specified.
#' The app allows the user to explore FSSA with simulated data, data that is provided on the server, or data that the user provides.
#' It allows the user to change parameters as they please, gives visual results of the methods, and also allows the user to compare FSSA results to other
#' spectrum analysis methods such as multivariate singular spectrum analysis. The tool is easy to use and can act as a nice starting point for a user that wishes to
#' perform FSSA as a part of their data analysis.
#'
#'
#'@seealso
#'  \code{\link{fssa}}, \code{\link{freconstruct}},
#'  \code{\link{fwcor}}, \code{\link{wplot}}, \code{\link{fts}}, \code{\link{plot.fts}}, \code{\link{plot.fssa}},
#'  \code{\link{cor.fts}}, \code{\link{launchApp}}
#'
#'
#'
#'@references
#'   Haghbin H., Najibi, S.M., Mahmoudvand R., Trinka J., Maadooliat M. (2019).
#'   Functional singular spectrum Analysis. Manuscript submitted for publication.
#'
#'
#'
#' @docType package
#' @name Rfssa
NULL
