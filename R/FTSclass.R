#' Functional Time Series Class
#'
#' This function is used to create functional time series (\code{\link{fts}}) objects from lists of discretely sampled data, basis specifications, and grid elements which provide the domain that each variable is observed over. Each variable is assumed to be observed over a regular and equidistant grid. In addition, each variable in the fts is assumed to be observed over a one or two-dimensional domain.
#' @param X  A list of length p where p is a positive integer specifying the number of variables observed in the fts. Each entry in the list should be a matrix or an array.
#' @param B A list of length p. Each entry in the list should be either a matrix specifying the basis for each variable or each list entry should be a list specifying the number of basis elements and desired basis type to be used in the smoothing process.
#' @param time A character vector where each entry specifies the time an observation is made.
#' @param method determine the `X` matrix type as "coefs" and "data".
#' @param argval A list of length p. Each entry in the list should either be a numeric or a list of numeric elements depending on the dimension of the domain the variable is observed over.
#' @importFrom fda eval.basis
#' @importFrom methods is new
#' @seealso \code{\link{fssa}}
#' @note Refer to \code{\link{fssa}} for examples on how to run this function.
#' @export
FTS <- function(X, B, argval = NULL, method = "data", start = 1, end = NULL ) { # Constructor function for the FTS class
  # Check if X is a matrix, and if so, convert it to a list
  if (is.array(X)) X <- list(X)
  if (is.basis(B) | is.array(B)) B <- list(B)
  if (is.numeric(argval)) argval <- list(argval)
  if (!is.list(X)) stop("The data input `X` must be a `matrix`,`array` or a list.")
  if (!is.list(B)) stop("The basis input `B` must be a `basisfd` object or a list of `basisfd` objects or empirical basis matrix.")
  # Check if all elements of X are matrices
  if (!all(sapply(X, is.array))) stop("All elements of the list X must be numeric `matrix` objects.")
  p <- length(X)
  dimSupp <- arg <- B_mat <- Coefs <- list()
  n_def <- 100
  for (j in 1L:p) {
    # determine the dimension support of the variable j
    dimSupp[[j]] <- ifelse(!is.basis(B[[j]]) && !is.array(B[[j]]), length(B[[j]]), 1)
    # Generating basis matrices=========================================
    # Generating a fd basis for variables whose domain is one-dimensional using a supplied list.
    if (dimSupp[[j]] == 1 && is.basis(B[[j]])) {
      # setting up the grids:
      if (is.null(argval)) {
        minval <- B[[j]]$rangeval[1]
        maxval <- B[[j]]$rangeval[2]
        if (method == "data") {
          arg[[j]] <- seq(from = minval, to = maxval, length.out = nrow(X[[j]]))
        } else { # method == "coefs" and NULL grids:
          arg[[j]] <- seq(from = minval, to = maxval, length.out = n_def)
        }
      } else {
        arg[[j]] <- argval[[j]]
      }
      B_mat[[j]] <- eval.basis(arg[[j]], basisobj = B[[j]])
    } else if (dimSupp[[j]] > 1) {
      if (all(sapply(B[[j]], is.basis))) {
        if (is.null(argval)) { # 2d fd basis with NULL grids
          minval1 <- B[[j]][[1]]$rangeval[1]
          maxval1 <- B[[j]][[1]]$rangeval[2]
          minval2 <- B[[j]][[2]]$rangeval[1]
          maxval2 <- B[[j]][[2]]$rangeval[2]
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
        u <- seq(1, nrow(B[[j]][[1]]))
        v <- seq(1, nrow(B[[j]][[2]]))
      } else { # Empirical 2d basis and not NULL grids
        u <- argval[[j]][[1]]
        v <- argval[[j]][[2]]
      }
      # Create the grid using u and v
      arg[[j]] <- cbind(rep(u, each = length(v)), rep(v, length(u)))
      b_1 <- eval.basis(evalarg = u, basisobj = B[[j]][[1]])
      b_2 <- eval.basis(evalarg = v, basisobj = B[[j]][[2]])
      B_mat[[j]] <- kronecker(b_1, b_2)

      M_x <- length(u)
      M_y <- length(v)
      # # Reshape FTS of Images from matricies to vectors
      if (is.matrix(X[[j]])) X[[j]] <- array(X[[j]], dim = c(M_x, M_y, 1))
      X[[j]] <- matrix(aperm(X[[j]], c(2, 1, 3)), nrow = M_x * M_y)
    }
    if (method == "data") {
      # Estimate the coefficients of each FTS variables.=========================================
      Coefs[[j]] <- solve(t(B_mat[[j]]) %*% B_mat[[j]]) %*% t(B_mat[[j]]) %*% X[[j]]
    } else { # method == "coefs"
      Coefs[[j]] <- X[[j]]
    }
  }
  # Creating the time object ========================================
  N <- tail(dim(X[[1]]), 1)
  if (is.null(end)) end <- start+N-1
  time <- seq(from = start, to = end, length.out = N)

  # Create and return an instance of the FTS class=========================================
  out <- list(Class = "FTS", N = N, dimSupp = dimSupp, time = time, Coefs = Coefs, Basis = B_mat, argval = arg)
  class(out) <- "FTS"
  return(out)
}

