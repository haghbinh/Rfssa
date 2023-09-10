#' Convert Object to a funts
#'
#' This function allows you to convert various types of objects into a functional time series (\code{\link{funts}}) object.
#'
#' @param obj The object to be converted. It can be an object of class \code{\link{fd}} (functional data) of the package `fda`, \code{\link{fts}} (functional time series) of the package `rainow`types.
#' @param basis An optional argument specifying the basis to be used for the resulting \code{\link{funts}} object when converting from \code{\link{fts}} objects. If not provided, a B-spline basis will be created by default.
#'
#' @return An object of class \code{\link{funts}}.
#'
#' @importFrom fda is.fd create.bspline.basis
#' @importFrom ftsa is.fts
#'
#' @examples
#' require(rainbow)
#' fts_obj=fts(x = 15:49,
#'             y = Australiasmoothfertility$y[,1:60],
#'             xname = "Age", yname = "Fertility rate")
#'
#' x_funts1 <- as.funts(fts_obj)
#' plot(x_funts1,
#'      main = "Australians Fertility",
#'      ylab = "Fertility rate",
#'      xlab = "Age")
#'
#' require(fda)
#' bs <- create.bspline.basis(rangeval = c(15,49),nbasis = 13)
#' fd_obj <- smooth.basis(argvals = fts_obj$x,fts_obj$y,bs)$fd
#'
#' x_funts <- as.funts(fd_obj)
#' plotly_funts(x_funts,
#'              main = "Australians Fertility",
#'              ylab = "Fertility rate",
#'              xlab = "Age")
#'
#' @seealso \code{\link{funts}}, \code{\link{create.bspline.basis}}, \code{\link{create.fts}}
#'
#' @note Only objects of class \code{\link{fd}} (functional data) and \code{\link{fts}} (functional time series) can be converted to a \code{\link{funts}} object. Other types will result in an error.
#'
#' @export
as.funts <- function(obj, basis = NULL) {
  if (is.fd(obj)) {
    x_funts <- funts(X = obj$coefs, basisobj = obj$basis, method = "coefs")
  } else if (is.fts(obj)) {
    argval <- obj$x
    m1 <- min(argval)
    m2 <- max(argval)
    X <- obj$y
    N <- ncol(X)
    time_st <- obj$time[1]
    time_en <- obj$time[N]
    d <- floor(min(dim(X)) / 2)
    if (is.null(basis)) basisobj <- create.bspline.basis(c(m1, m2), nbasis = d)
    x_funts <- funts(X = X, basisobj = basisobj, argval = argval, start = time_st, end = time_en)
  } else {
    stop("This object can not be converted to a funts object.")
  }
  return(x_funts)
}
