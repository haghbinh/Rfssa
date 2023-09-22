#' Rfssa: A Package for Functional Singular Spectrum Analysis and Related Methods.
#'
#' @description
#' The Rfssa package provides essential functions for implementing functional singular spectrum analysis (FSSA) methods. It supports the analysis of univariate and multivariate functional time series (FTS). FSSA, both univariate and multivariate, is a non-parametric approach for decomposing and reconstructing FTS. This package allows you to work with FTS variables observed over one or two-dimensional domains.
#'
#' @details
#' To use this package, start by creating a \code{\link{funts}} object. You can
#' define it by providing raw data, basis specifications, and grid specifications.
#'  The FTS object can be univariate or multivariate, and the variables can be observed
#'  over one-dimensional curves or two-dimensional images. The package includes
#'  input validity checks to guide you.
#' Use the \code{\link{plot}} method to visualize the \code{\link{funts}} object.
#'  It offers various plotting options for one-dimensional domain variables and
#'  animations for two-dimensional domains.
#' Next, apply the FSSA routine (\code{\link{fssa}}) to the \code{\link{funts}}
#' object with a chosen lag parameter to obtain the decomposition. The decomposition
#'  function utilizes `RSpectra` and `RcppEigen` R packages, along with the `Eigen` C++ package,
#'   for efficient processing.
#' Visualize the decomposition results using the \code{\link{plot.fssa}} method to
#'  help choose the grouping of eigentriples for reconstruction (\code{\link{freconstruct}})
#'   or forecasting (\code{\link{fforecast}}). The \code{\link{freconstruct}}
#'   function reconstructs a list of \code{\link{funts}} objects specified by the
#'   grouping, while \code{\link{fforecast}} provides predictions of the signals
#'   specified by the grouping. Calculate bootstrapped prediction intervals for forecasts
#'   using the \code{\link{fpredinterval}} function.
#' When performing forecasting, typically specify one group that captures the assumed
#' deterministic, extracted signal within the FTS, while excluding all other modes of variation.
#'  Currently, forecasting supports FTS with one-dimensional domains, with two-dimensional
#'   domain forecasting planned for future updates.
#' This version of the `Rfssa` R package introduces the \code{\link{fpredinterval}} function,
#' which calculates prediction intervals for FSSA-based forecasts using a bootstrap
#' approach for the residuals.
#'
#' @seealso
#'  \code{\link{fssa}}, \code{\link{freconstruct}}, \code{\link{fforecast}}
#'  \code{\link{funts}}, \code{\link{launchApp}}
#'
#'
#' @references
#'   Haghbin, H., Najibi, S. M., Mahmoudvand, R., Trinka, J., Maadooliat, M. (2021). Functional singular spectrum analysis. Stat, 10(1), e330.
#'
#'   Trinka J. (2021) Functional Singular Spectrum Analysis: Nonparametric Decomposition and Forecasting Approaches for Functional Time Series [Doctoral dissertation, Marquette University]. ProQuest Dissertations Publishing.
#'
#'   Trinka, J., Haghbin, H.,  Maadooliat, M. (2022). Multivariate Functional Singular Spectrum Analysis: A Nonparametric Approach for Analyzing Multivariate Functional Time Series. In Innovations in Multivariate Statistical Modeling: Navigating Theoretical and Multidisciplinary Domains (pp. 187-221). Cham: Springer International Publishing.
#'
#'   Trinka, J., Haghbin, H., Shang, H., Maadooliat, M. (2023). Functional Time Series Forecasting: Functional Singular Spectrum Analysis Approaches. Stat, e621.
#'

#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib Rfssa, .registration = TRUE
#' @importFrom fda is.fd create.bspline.basis eval.basis is.basis pca.fd eval.fd
#' @importFrom rainbow fts
#' @importFrom ftsa quantile.fts
#' @importFrom RSpectra eigs
#' @importFrom lattice xyplot levelplot
#' @importFrom plotly renderPlotly plotlyOutput plot_ly add_lines layout subplot add_surface hide_colorbar ggplotly
#' @importFrom graphics abline axis par plot points polygon
#' @importFrom stats fft integrate rnorm sd ts.plot density approxfun
#' @importFrom utils head read.table
#' @importFrom markdown mark
#' @importFrom Rssa ssa reconstruct
#' @importFrom ggplot2 ggplot aes_string unit geom_tile scale_fill_distiller xlab ylab labs ggtitle scale_y_continuous scale_x_continuous waiver theme element_line element_text element_blank
#' @importFrom tibble as_tibble
#' @importFrom dplyr "%>%" group_by
## usethis namespace: end
NULL
