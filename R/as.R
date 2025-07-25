#' Convert Object to a funts
#'
#' This function allows you to convert various types of objects into a functional time series (\code{\link{funts}}) object.
#'
#' @param obj the object to be converted. It can be an object of class \code{\link[fda]{fd}} (functional data) of the package `fda`, \code{\link[rainbow]{fts}} (functional time series) of the package `rainbow` types.
#' @param basis an optional argument specifying the basis to be used for the resulting \code{\link{funts}} object when converting from \code{\link[rainbow]{fts}} objects. If not provided, a B-spline basis will be created by default.
#'
#' @return An object of class \code{\link{funts}}.
#'
#'
#' @examples
#' require(rainbow)
#' class(Australiasmoothfertility)
#' x_funts1 <- as.funts(Australiasmoothfertility)
#' plot(x_funts1, main = "Australians Fertility")
#'
#' require(fda)
#' bs <- create.bspline.basis(rangeval = c(15, 49), nbasis = 13)
#' fd_obj <- smooth.basis(argvals = Australiasmoothfertility$x, Australiasmoothfertility$y, bs)$fd
#'
#' x_funts <- as.funts(fd_obj)
#' plotly_funts(x_funts,
#'   main = "Australians Fertility",
#'   ylab = "Fertility rate",
#'   xlab = "Age"
#' )
#' @seealso \code{\link{funts}}, \code{\link[fda]{create.bspline.basis}}
#'
#' @note Only objects of class \code{\link[fda]{fd}} (functional data) and \code{\link[rainbow]{fts}} (functional time series) can be converted to a \code{\link{funts}} object. Other types will result in an error.
#'
#' @export
as.funts <- function(obj, basis = NULL) {
  if (is.fd(obj)) {
    x_funts <- funts(X = obj$coefs, basisobj = obj$basis, method = "coefs")
  } else if (inherits(obj, c("fds", "fts"))) {
    argval <- obj$x
    m1 <- min(argval)
    m2 <- max(argval)
    X <- obj$y
    N <- ncol(X)
    time_st <- 1
    time_en <- NULL
    if (inherits(obj, "fts")) {
      time_st <- obj$time[1]
      time_en <- obj$time[N]
    }
    vnames <- obj$yname
    dnames <- list(obj$xname)
    d <- floor(min(dim(X)) / 2)
    if (is.null(basis)) basisobj <- create.bspline.basis(c(m1, m2), nbasis = d)
    x_funts <- funts(X = X, basisobj = basisobj, argval = argval, start = time_st,
                     end = time_en, vnames = vnames, dnames = dnames)
  } else {
    stop("This object can not be converted to a funts object.")
  }
  return(x_funts)
}
