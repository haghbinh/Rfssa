#' Rfssa: A Package for Functional Singular Spectrum Analysis and Related Methods.
#'
#' The Rfssa package provides essential functions for implementing functional singular spectrum analysis (FSSA) methods. It supports the analysis of univariate and multivariate functional time series (FTS). FSSA, both univariate and multivariate, is a non-parametric approach for decomposing and reconstructing FTS. This package allows you to work with FTS variables observed over one or two-dimensional domains.
#'
#' @details
#'
#' To use this package, start by creating a \code{\link{funts}} object. You can define it by providing raw data, basis specifications, and grid specifications. The FTS object can be univariate or multivariate, and the variables can be observed over one-dimensional curves or two-dimensional images. The package includes input validity checks to guide you.
#'
#' Use the \code{\link{plot.funts}} method to visualize the \code{\link{funts}} object. It offers various plotting options for one-dimensional domain variables and animations for two-dimensional domains.
#'
#' Next, apply the FSSA routine (\code{\link{fssa}}) to the \code{\link{funts}} object with a chosen lag parameter to obtain the decomposition. The decomposition function utilizes RSpectra and RcppEigen R packages, along with the Eigen C++ package, for efficient processing.
#'
#' Visualize the decomposition results using the \code{\link{plot.fssa}} method to help choose the grouping of eigentriples for reconstruction (\code{\link{freconstruct}}) or forecasting (\code{\link{fforecast}}). The \code{\link{freconstruct}} function reconstructs a list of \code{\link{funts}} objects specified by the grouping, while \code{\link{fforecast}} provides predictions of the signals specified by the grouping. Calculate bootstrapped prediction intervals for forecasts using the \code{\link{fpredinterval}} function.
#'
#' When performing forecasting, typically specify one group that captures the assumed deterministic, extracted signal within the FTS, while excluding all other modes of variation. Currently, forecasting supports FTS with one-dimensional domains, with two-dimensional domain forecasting planned for future updates.
#'
#' This version of the Rfssa R package introduces the \code{\link{fpredinterval}} function, which calculates prediction intervals for FSSA-based forecasts using a bootstrap approach for the residuals.
#'
#' @seealso
#'  \code{\link{fssa}}, \code{\link{freconstruct}}, \code{\link{fforecast}}
#'  \code{\link{fwcor}}, \code{\link{wplot}}, \code{\link{funts}}, \code{\link{plot.funts}}, \code{\link{plot.fssa}},
#'  \code{\link{launchApp}}
#'
#'
#' @references
#'   Haghbin, H., Najibi, S. M., Mahmoudvand, R., Trinka, J., Maadooliat, M. (2021). Functional singular spectrum analysis. Stat, 10(1), e330.
#'
#'   Trinka J., Haghbin H., Maadooliat M. (Accepted) Multivariate Functional Singular Spectrum Analysis: A Nonparametric Approach for Analyzing Multivariate Functional Time Series. In: Bekker A., Ferreira, J., Arashi M., Chen D. (eds) Innovations in Multivariate Statistical Modeling: Navigating Theoretical and Multidisciplinary Domains. Emerging Topics in Statistics and Biostatistics. Springer, Cham.
#'
#'   Trinka J. (2021) Functional Singular Spectrum Analysis: Nonparametric Decomposition and Forecasting Approaches for Functional Time Series [Doctoral dissertation, Marquette University]. ProQuest Dissertations Publishing.
#'
#'   Trinka, J., Haghbin, H.,  Maadooliat, M. (2022). Multivariate Functional Singular Spectrum Analysis: A Nonparametric Approach for Analyzing Multivariate Functional Time Series. In Innovations in Multivariate Statistical Modeling: Navigating Theoretical and Multidisciplinary Domains (pp. 187-221). Cham: Springer International Publishing.
#'
#'   Trinka, J., Haghbin, H., Shang, H., Maadooliat, M. (2023). Functional Time Series Forecasting: Functional Singular Spectrum Analysis Approaches. Stat, DOI: 10.1002/sta4.621
#'
#' @docType package
#' @name Rfssa
NULL
