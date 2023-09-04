#' Length of Functional Time Series
#'
#' Returns the length of a "funts" object.
#' @param obj
#' @export
length.funts <- function(obj) {
  return(obj$N)
}
#=======================================================================

#' Custom Print Method for Functional Time Series (funts) Objects
#'
#' This function provides a custom print method for objects of the Functional Time Series (funts) class.
#' It displays information about the funts object, including its length (N), the number of variables,
#' and the structure of the object.
#'
#' @param obj An object of class "funts" to be printed.
#'
#' @export
print.funts <- function(obj) {
  cat('\nFunctional time series (funts) object:')
  cat("\nlength N: ", obj$N)
  cat("\nnumber of variables: ", length(obj$dimSupp),'\n')
  cat(str(obj),'\n')
}
#=======================================================================

#' Addition of Functional Time Series
#'
#' A method for functional time series (\code{\link{funts}}) addition and funts-scalar addition. Note that if the funts is multivariate
#' then a vector of numerics may be provided allowing for addition of different scalars to different variables. For example, multivariate funts-numeric addition
#' follows the form of \code{Y+c(1,2)} if \code{Y} is a bivariate FTS.
#' @return An object of class \code{\link{funts}}.
#'
#' @param Y1 An object of class \code{\link{funts}} or numeric.
#' @param Y2 An object of class \code{\link{funts}} or numeric.

#' @seealso \code{\link{funts}}
#' @export
"+.funts" <- function(obj1, obj2 = NULL) {
  if (is.null(obj2)) {
    return(obj1)
  }
  if (xor(is.funts(obj1), is.funts(obj2))) {
    if (!xor(is.double(obj1), is.double(obj2))) stop("At least one object must be an funts, and the other one can be a scalar")
    if (is.double(obj1)) {
      temp <- obj1
      obj1 <- obj2
      obj2 <- temp
    }
    coef <- lapply(obj1$coefs, function(x) x * obj2)
  } else { # if (is.funts(Y1) && is.funts(Y2))
    if (length(obj1$coefs) != length(obj2$coefs)) stop("Error: Incompatible Data Lengths. Functional time series require data of the same length.")
    if (!all.equal(obj1$time, obj2$time)) stop("Error: Incompatible Data Time index. Functional time series require data of the same Time index.")
    if (!all.equal(obj1$basis, obj2$basis)) stop("Error: Incompatible Data basis. Functional time series require the same basis functions.")
    if (!all.equal(obj1$dimSupp, obj2$dimSupp)) stop("Error: Incompatible Data dimension supports. Functional time series require the same dimSupps.")
    coef <- obj1$coefs + obj2$coefs
  }
  obj1$coefs <- coef
  return(obj1)
}
#=======================================================================

#' Scalar Multiplication of a Functional Time Series (\code{\link{funts}}) Object
#'
#' Performs scalar multiplication of a Functional Time Series (funts) object by either another funts object or a scalar value.
#'
#' @param obj1 An object of class \code{\link{funts}} or a scalar value.
#' @param obj2 An object of class \code{\link{funts}} or a scalar value.
#' @return An object of class \code{\link{funts}} representing the result of scalar multiplication.
#'
#' @details This function allows you to multiply a Functional Time Series (funts) object by either another funts object or a scalar value.
#' The operation is element-wise when both objects are funts, and scalar multiplication when one object is a funts and the other is a scalar.
#'
#' @seealso \code{\link{funts}}
#'
#' @export
"*.funts" <- function(obj1, obj2) {
  if (xor(is.funts(obj1), is.funts(obj2))) {
    if (xor(is.double(obj1), is.double(obj2))) {
      if (is.double(obj1)) {
        temp <- obj1
        obj1 <- obj2
        obj2 <- temp
      }
    }
    obj1$coefs <- lapply(obj1$coefs, function(x) x * obj2)
  } else {
    stop("One object must be an funts, and the other one a scalar")
  }
  return(obj1)
}
#=======================================================================
#'
#' Subtract two `funts` objects or a `funts` object and a scalar
#'
#' This function allows you to subtract two `funts` (functional time series) objects or
#' subtract a `funts` object from a scalar. It returns the difference between the two
#' objects.
#'
#' @param obj1 An object of class `funts`. This represents the first `funts` object.
#' @param obj2 An object of class `funts` or a scalar.
#'
#' @seealso \code{\link{funts}}
#'
#' @export
"-.funts" <- function(obj1, obj2 = NULL) {
  # Check if at least one object must be an `funts`, and the other one can be a scalar
  if (!xor(is.funts(obj1), is.funts(obj2)) && !is.null(obj2)) {
    stop("At least one object must be a `funts`, and the other one can be a scalar")
  }
  if (is.null(obj2)) {
    # If obj2 is NULL, subtract obj1 from zero
    return((-1) * obj1)
  } else {
    # Subtract obj2 from obj1 element-wise
    return(obj1 + (-1) * obj2)
  }
}
#=======================================================================
#'



