#' Length of Functional Time Series
#'
#' Returns the length of a "funts" object.
#' @param obj
#' @export
length.funts <- function(obj){
  return(obj$N)
}


#'
#' Addition of Functional Time Series
#'
#' A method for functional time series (\code{\link{funts}}) addition and funts-scalar addition. Note that if the funts is multivariate
#' then a vector of numerics may be provided allowing for addition of different scalars to different variables. For example, multivariate funts-numeric addition
#' follows the form of \code{Y+c(1,2)} if \code{Y} is a bivariate FTS.
#' @return An object of class \code{\link{funts}}.
#'
#' @param Y1 An object of class \code{\link{funts}} or numeric.
#' @param Y2 An object of class \code{\link{funts}} or numeric.
#' @examples
#' \dontrun{
#' require(Rfssa)
#' load_github_data("https://github.com/haghbinh/Rfssa/blob/master/data/Callcenter.RData")
#' D <- matrix(sqrt(Callcenter$calls), nrow = 240)
#' u <- seq(0, 1, length.out = 240) # Define domain of functional data
#' d <- 22 # number of basis elements
#' Y <- Rfssa::funts(list(D), list(list(d, "bspline")), list(u))
#' plot(Y)
#' Yplus <- Y + Y # add the functional time series to itself
#' plot(Yplus)
#' Yplus2 <- Y + 2 # add 2 to every term in the functional time series
#' plot(Yplus2)
#' }
#' @seealso \code{\link{funts}}
#' @useDynLib Rfssa
#' @export
"+.funts" <- function(Y1, Y2) {
  if (is.funts(Y1) && is.funts(Y2)) {
    if (length(Y1) != length(Y2)) stop("Error: Incompatible Data Lengths. Functional time series require data of the same length.")
    if (!all.equal(Y1$time, Y2$time)) stop("Error: Incompatible Data Time index. Functional time series require data of the same Time index.")
    if (!all.equal(Y1$basisobj, Y2$basisobj)) stop("Error: Incompatible Data basis. Functional time series require of the same basis functions.")
  }
  time <- Y1$time
  p <- length(Y1$dimSupp)
  grid <- Y1$argval
  basis <- Y1$Basis
  eval_funts <- list()
  new_grid <- list()

}
