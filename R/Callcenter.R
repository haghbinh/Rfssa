#' Callcenter Dataset: Number of Calls for a Bank
#'
#' This dataset represents a small call center for an anonymous bank (Brown et al., 2005).
#' It provides detailed information about the exact times of calls that were connected
#' to the center throughout the year 1999, from January 1 to December 31.
#'
#' The data have been converted into a functional time series using a B-spline basis
#' system with 22 basis functions. The resulting dataset is stored as a functional time
#' series object of class `funts`. You can load the raw data using the function
#' \code{\link{loadCallcenterData}}. See \code{\link{funts}} for more details.
#'
#' @name Callcenter
#' @format A functional time series object of class `funts` with the following fields:
#' \describe{
#'   \item{time}{the time index indicating when the calls occurred.}
#'   \item{coefs}{the coefficients corresponding to the B-spline basis functions.}
#'   \item{basisobj}{the basis functions used for the functional representation.}
#'   \item{dimSupp}{the dimension support of the functional data.}
#' }
#'
#' @seealso \code{\link{loadCallcenterData}}, \code{\link{funts}}
#'
#' @references
#' \enumerate{
#' \item
#' Brown, L., Gans, N., Mandelbaum, A., Sakov, A., Shen, H., Zeltyn, S., & Zhao, L. (2005).
#' Statistical analysis of a telephone call center: A queueing-science perspective.
#' \emph{Journal of the American Statistical Association}, \strong{100}(469), 36-50.
#' }
#'
#' @examples
#' \dontrun{
#' # Load the Callcenter dataset
#' data("Callcenter")
#' }
NULL
