#'
#' Length of Functional Time Series
#'
#' Returns the length of a "funts" object.
#' @param x an object of class "funts" .
#'
#' @examples
#' data("Callcenter")
#' length(Callcenter)
#'
#' @export
length.funts <- function(x) {
  return(x$N)
}
# =======================================================================
#' Custom Print Method for Functional Time Series (funts) Objects
#'
#' This custom print method is designed for objects of the Functional Time Series (funts) class.
#' It provides a summary of the funts object, including its length (N), the number of variables,
#' and its structure.
#'
#' @param x an object of class "funts" to be printed.
#' @param ...	further arguments passed to or from other methods.
#' @examples
#' data("Callcenter")
#' print(Callcenter)
#'
#' @export
print.funts <- function(x, ...) {
  cat("\nFunctional time series (funts) object:")
  cat("\nNumber of variables: ", length(x$dimSupp))
  cat("\nLenght: ", x$N)
  cat("\nStart: ", x$time[1])
  cat("\nEnd: ", x$time[x$N])
  cat("\nTime: ")
  cat(str(x$time))
}
# =======================================================================
#'
#' Addition of Functional Time Series
#'
#' A method for functional time series (\code{\link{funts}}) addition and funts-scalar addition.
#' @return an object of class \code{\link{funts}}.
#'
#' @param obj1 an object of class \code{\link{funts}} or numeric.
#' @param obj2 an object of class \code{\link{funts}} or numeric.
#'
#'
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
    coef <- list()
    for (i in 1:length(obj1$coefs)) coef[[i]] <- obj1$coefs[[i]] + obj2$coefs[[i]]
  }
  obj1$coefs <- coef
  return(obj1)
}
# =======================================================================
#' Scalar Multiplication of Functional Time Series (funts)
#'
#' Perform scalar multiplication of a Functional Time Series (funts) object by either another funts object or a scalar value.
#'
#' @param obj1 an object of class \code{\link{funts}} or a scalar value.
#' @param obj2 an object of class \code{\link{funts}} or a scalar value.
#' @return An object of class \code{\link{funts}} representing the result of scalar multiplication.
#'
#' @details This function allows element-wise multiplication of a Functional Time Series (funts) object by another funts object or scalar value.
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
    stop("One object must be an funts, and the other one a scalar.")
  }
  return(obj1)
}
# =======================================================================
#'
#' Subtract `funts` Objects or a `funts` Object and a Scalar
#'
#' This function allows subtraction between two `funts` (functional time series) objects or
#' between a `funts` object and a scalar. It returns the resulting difference.
#'
#' @param obj1 an object of class `funts`, representing the first `funts` object.
#' @param obj2 an object of class `funts` or a scalar.
#'
#' @seealso \code{\link{funts}}
#'
#' @examples
#' data("Callcenter")
#' y <- Callcenter
#' print(1 - y)
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
# =======================================================================
#'
#' Indexing into Functional Time Series
#'
#' An indexing method for functional time series (\code{\link{funts}}) objects.
#'
#' This function allows you to extract specific subsets of a functional time series
#' based on the provided indices. You can specify which subsets you want to extract
#' from the functional time series.
#'
#' @param obj an object of class \code{\link{funts}}.
#' @param i an index or indices specifying the subsets of times to extract.
#' @param j an index or indices specifying the subsets of variables to extract.
#' @return an `funts` object containing the specified subsets
#' @seealso \code{\link{funts}}
#' @export
"[.funts" <- function(obj, i = NULL, j = NULL) {
  p <- length(obj$dimSupp)
  if (is.null(i)) i <- 1:obj$N
  if (is.null(j)) j <- 1:p
  if (max(i) > obj$N | min(i) < 1) stop(" subscript i out of bounds")
  if (max(j) > p | min(j) < 1) stop(" subscript j out of bounds")
  out <- list(N = length(i), dimSupp = list(), time = obj$time[i], coefs = list(), basis = list(), B_mat = list(), argval = list())
  count <- 0
  for (k in j) {
    count <- count + 1
    out$dimSupp[[count]] <- obj$dimSupp[[k]]
    out$basis[[count]] <- obj$basis[[k]]
    out$B_mat[[count]] <- obj$B_mat[[k]]
    out$argval[[count]] <- obj$argval[[k]]
    out$coefs[[count]] <- as.matrix(obj$coefs[[k]])[, i]
  }
  class(out) <- "funts"
  return(out)
}
# =======================================================================
#'
#' Evaluate a Functional Time Series (funts) Object on a Given Grid
#'
#' This function allows you to evaluate a Functional Time Series (funts) object
#' on a specified grid of argument values. The result is a list of matrices, each
#' matrix corresponding to one dimension of the functional data.
#'
#' @param argvals a list or numeric vector specifying the grid points at which
#'   to evaluate the functional time series. For multivariate functional data,
#'   provide a list of grids corresponding to each dimension.
#' @param obj an object of class \code{\link{funts}} to be evaluated.
#'
#' @return A list of matrices, where each matrix represents the evaluated values
#'   of the functional data on the specified grid.
#'
#' @details
#' The \code{argvals} argument can be a list of grids for multivariate functional data.
#' The function handles both functional basis and empirical basis cases for evaluation.
#' For empirical basis with irregular grids, a warning is issued as this feature
#' is under development.
#'
#' @seealso \code{\link{funts}}
#'
#' @examples
#' data("Montana")
#' y <- Montana
#' u <- seq(0, 23, len = 4)
#' v <- seq(1, 33, len = 3)
#' grid <- list(u, list(v, v))
#' eval.funts(grid, y)
#' @export
eval.funts <- function(argvals, obj) {
  if (!xor(is.numeric(argvals), is.list(argvals))) stop("Error: Incompatible grid points. It must be a list or numeric object.")
  if (!is.funts(obj)) stop("The obj argument must have class of funts.")
  dimSupp <- obj$dimSupp
  p <- length(dimSupp)
  basis <- obj$basis
  result <- list()
  for (j in 1:p) {
    if (is.numeric(argvals) || is.array(argvals) || (p == 1 && length(argvals) > 1)) argvals <- list(argvals)
    if (dimSupp[[j]] == 1) {
      if (is.basis(basis[[j]])) {
        B <- eval.basis(evalarg = argvals[[j]], basisobj = basis[[j]])
      } else { # if(is.matrix(basis[[j]])))  #Empirical basis
        B <- eval.empb(evalarg = argvals[[j]], basisobj = basis[[j]])
      }
      result[[j]] <- B %*% obj$coefs[[j]]
    } else { # dimSupp[[j]]==2
      u <- argvals[[j]][[1]]
      v <- argvals[[j]][[2]]
      if (all(sapply(basis[[j]], is.basis))) { # fd basis
        b_1 <- eval.basis(evalarg = u, basisobj = basis[[j]][[1]])
        b_2 <- eval.basis(evalarg = v, basisobj = basis[[j]][[2]])
      } else { # Empirical basis (list)
        if (is.list(basis[[j]])) {
          b_1 <- eval.empb(evalarg = u, basisobj = basis[[j]][[1]])
          b_2 <- eval.empb(evalarg = v, basisobj = basis[[j]][[2]])
        } else { # Empirical basis (Kronecker product or SVD)
          warning("Warning: This feature is under development for the basis defined over irregular grids. Results may not be accurate.")
          b_1 <- basis[[j]][, ]
          b_2 <- 1
        }
      }
      B <- kronecker(b_1, b_2)
      result[[j]] <- aperm(array(B %*% obj$coefs[[j]], dim = c(length(v), length(u), obj$N)), c(2, 1, 3))
    }
  }
  return(result)
}
# =======================================================================
#' Plot Functional Time Series (funts) Data
#'
#' Create visualizations of Functional Time Series (funts) data, supporting both one-dimensional and two-dimensional domains.
#'
#' @param x An object of class \code{funts}.
#' @param npts Number of grid points for the plots.
#' @param obs Observation number (for two-dimensional domains).
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param main Main title for the plot.
#' @param type Type of plot ("l" for line, "p" for points, etc.).
#' @param lty Line type (1 for solid, 2 for dashed, etc.).
#' @param ... Additional graphical parameters passed to plotting functions.
#'
#' @details
#' This function enables the creation of visualizations for Functional Time Series (funts) data. It supports both one-dimensional and two-dimensional domains.
#'
#' For one-dimensional domains, line plots are used, while for two-dimensional domains, image plots are employed.
#'
#'
#' @importFrom graphics par image axis
#'
#' @examples
#' # Example with one-dimensional domain
#' data("Callcenter")
#' plot(Callcenter,
#'   lwd = 2, col = "deepskyblue4",
#'   main = "Call Center Data",
#'   xlab = "Time (6 minutes aggregated)",
#'   ylab = "Sqrt of Call Numbers"
#' )
#'
#' # Example with two-dimensional domain
#' data("Montana")
#' plot(Montana,
#'   obs = 2,
#'   main = c("Temperature Curves", "NDVI Images,"),
#'   xlab = c("Time", "Longitude"),
#'   ylab = c("Normalized Temperature (\u00B0C)", "Latitude")
#' )
#'
#' @seealso \code{\link{funts}}, \code{\link{Callcenter}}, \code{\link{Montana}}
#' @importFrom graphics matplot
#'
#' @export
plot.funts <- function(x, npts = 100, obs = 1, xlab = NULL, ylab = NULL, main = NULL, type = "l", lty = 1, ...) {
  dimSupp <- x$dimSupp
  p <- length(dimSupp)
  if (is.null(xlab)) {
    for (i in 1:p){
      xlab[i] <- x$dnames[[i]][1]
    }
  }
  if (is.null(ylab)) {
    for (i in 1:p){
      if (dimSupp[[i]] == 1) {
        ylab[i] <- x$vnames[i]
      } else { # dimSupp[[i]] == 2
        ylab[i] <- x$dnames[[i]][2]
      }
    }
  }
  N <- x$N
  time <- x$time
  par(mfrow = c(p, 1))
  if (is.null(ylab)) ylab <- paste("Variable ", 1:p)
  if (is.null(xlab)) xlab <- rep("time", p)
  for (j in 1:p) {
    if (dimSupp[[j]] == 1) {
      supp <- matrix(range(x$argval[[j]]), nrow = 2)
      x_grids <- seq(supp[1, 1], supp[2, 1], len = npts)
      X <- eval.funts(x_grids, x[, j])[[1]]
      matplot(x_grids, X, type = type, lty = lty, xlab = xlab[j], ylab = ylab[j], main = main[j], ...)
    } else { # dim >1
      rangeval <- apply(x$argval[[j]], 2, range)
      supp <- matrix(c(rangeval[, 1], rangeval[, 2]), nrow = 2)
      x_grids <- seq(supp[1, 1], supp[2, 1], len = npts)
      y_grids <- seq(from = supp[1, 2], to = supp[2, 2], length.out = npts)
      grids2d <- list(x_grids, y_grids)
      X <- eval.funts(grids2d, x[, j])[[1]]
      image(X[, , obs], xlab = xlab[j], ylab = ylab[j], axes = FALSE, main = paste(main[j], " Observation:", obs), ...)
      axis(side = 1, at = seq(from = 0, to = 1, length.out = 10), labels = round(seq(supp[1, 1], supp[2, 1], len = 10), 1))
      axis(side = 2, at = seq(from = 0, to = 1, length.out = 10), labels = round(seq(supp[1, 2], supp[2, 2], len = 10), 1))
    }
  }
}
