#' Multivariate Functional Time Series (mfts) Class
#'
#' The `mmfts` class is designed to encapsulate functional time series objects, both univariate (FTS) and multivariate (MFTS).
#' It allows you to create `mfts` objects with specified basis systems. You can create an `mfts` object either by providing discrete
#' observed samples of curves and applying smoothing using a basis method or by supplying the set of basis coefficients.
#' `mfts` objects can handle two-dimensional functional data and different dimensions for multivariate FTS.
#'
#' @param X A matrix, three-dimensional array, or a list containing matrix or array objects. When `method="data"`, it represents
#' the set of curve values at discrete sampling points or argument values. When `method="coefs"`, `X` specifies the coefficients
#' corresponding to the basis system defined in 'basisobj'. If `X` is a list, it defines a multivariate FTS, with each element
#' being a matrix or three-dimensional array object. In matrix objects, rows correspond to argument values, and columns correspond
#' to the length of the FTS. In three-dimensional array objects, the first and second dimensions correspond to argument values,
#' and the third dimension to the length of the FTS.
#'
#' @param basisobj An object of class `basisfd`, a matrix of empirical basis, or a list of `basisfd` or empirical basis objects.
#' For empirical basis, rows correspond to basis functions, and columns correspond to grid points.
#'
#' @param argval A vector list of length `p`. Each entry in the list should either be a numeric value or a list of numeric
#' elements, depending on the dimension of the domain the variable is observed over.
#'
#' @param method Determines the type of the `X` matrix: "coefs" or "data."
#'
#' @param start The time of the first observation. It can be a single number or an object of classes `Date`, `POSIXct`, or `POSIXt`,
#' specifying a natural time unit.
#'
#' @param end The time of the last observation, specified in the same way as `start`.
#'
#' @importFrom fda eval.basis
#' @importFrom methods is new
#'
#' @seealso \code{\link{fssa}}
#'
#' @note Refer to \code{\link{fssa}} for examples on how to use this function.
#'
#' @export
mfts <- function(X, basisobj, argval = NULL, method = "data", start = 1, end = NULL ) { # Constructor function for the mfts class
  # Check if X is a matrix, and if so, convert it to a list
  if (is.array(X)) X <- list(X)
  if (is.basis(basisobj) | is.array(B)) basisobj <- list(basisobj)
  if (is.numeric(argval)) argval <- list(argval)
  if (!is.list(X)) stop("The data input `X` must be a `matrix`,`array` or a list.")
  if (!is.list(basisobj)) stop("The basis input `basisobj` must be a `basisfd` object or a list of `basisfd` objects or empirical basis matrix.")
  # Check if all elements of X are matrices
  if (!all(sapply(X, is.array))) stop("All elements of the list X must be numeric `matrix` objects.")
  p <- length(X)
  dimSupp <- arg <- B_mat <- Coefs <- list()
  n_def <- 100
  for (j in 1L:p) {
    # determine the dimension support of the variable j
    dimSupp[[j]] <- ifelse(!is.basis(basisobj[[j]]) && !is.array(basisobj[[j]]), length(basisobj[[j]]), 1)
    # Generating basis matrices=========================================
    # Generating a fd basis for variables whose domain is one-dimensional using a supplied list.
    if (dimSupp[[j]] == 1 && is.basis(B[[j]])) {
      # setting up the grids:
      if (is.null(argval)) {
        minval <- basisobj[[j]]$rangeval[1]
        maxval <- basisobj[[j]]$rangeval[2]
        if (method == "data") {
          arg[[j]] <- seq(from = minval, to = maxval, length.out = nrow(X[[j]]))
        } else { # method == "coefs" and NULL grids:
          arg[[j]] <- seq(from = minval, to = maxval, length.out = n_def)
        }
      } else {
        arg[[j]] <- argval[[j]]
      }
      B_mat[[j]] <- eval.basis(arg[[j]], basisobj = basisobj[[j]])
    } else if (dimSupp[[j]] > 1) {
      if (all(sapply(B[[j]], is.basis))) {
        if (is.null(argval)) { # 2d fd basis with NULL grids
          minval1 <- basisobj[[j]][[1]]$rangeval[1]
          maxval1 <- basisobj[[j]][[1]]$rangeval[2]
          minval2 <- basisobj[[j]][[2]]$rangeval[1]
          maxval2 <- basisobj[[j]][[2]]$rangeval[2]
          if (method == "data") {
            u <- seq(from = minval1, to = maxval1, length.out = dim(X[[j]])[2])
            v <- seq(from = minval2, to = maxval2, length.out = dim(X[[j]])[2])
          } else { # method == "coefs" and NULL grids:
            u <- seq(from = minval1, to = maxval1, length.out = n_def)
            v <- seq(from = minval2, to = maxval2, length.out = n_def)
          }
        } else { #  2d fd basis Not NULL grids:
          u <- argval[[j]][[1]]
          v <- argval[[j]][[2]]
        }
      } else if (is.null(argval)) { # Empirical 2d basis and NULL grids
        u <- seq(1, nrow(basisobj[[j]][[1]]))
        v <- seq(1, nrow(basisobj[[j]][[2]]))
      } else { # Empirical 2d basis and not NULL grids
        u <- argval[[j]][[1]]
        v <- argval[[j]][[2]]
      }
      # Create the grid using u and v
      arg[[j]] <- cbind(rep(u, each = length(v)), rep(v, length(u)))
      b_1 <- eval.basis(evalarg = u, basisobj = basisobj[[j]][[1]])
      b_2 <- eval.basis(evalarg = v, basisobj = basisobj[[j]][[2]])
      B_mat[[j]] <- kronecker(b_1, b_2)

      M_x <- length(u)
      M_y <- length(v)
      # # Reshape mfts of Images from matricies to vectors
      if (is.matrix(X[[j]])) X[[j]] <- array(X[[j]], dim = c(M_x, M_y, 1))
      X[[j]] <- matrix(aperm(X[[j]], c(2, 1, 3)), nrow = M_x * M_y)
    }
    if (method == "data") {
      # Estimate the coefficients of each mfts variables.=========================================
      Coefs[[j]] <- solve(t(B_mat[[j]]) %*% B_mat[[j]]) %*% t(B_mat[[j]]) %*% X[[j]]
    } else { # method == "coefs"
      Coefs[[j]] <- X[[j]]
    }
  }
  # Creating the time object ========================================
  N <- tail(dim(X[[1]]), 1)
  if (is.null(end)) end <- start+N-1
  time <- seq(from = start, to = end, length.out = N)

  # Create and return an instance of the mfts class=========================================
  out <- list(Class = "mfts", N = N, dimSupp = dimSupp, time = time, Coefs = Coefs, Basis = B_mat, argval = arg)
  class(out) <- "mfts"
  return(out)
}

