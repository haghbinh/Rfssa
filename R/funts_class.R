#' Functional Time Series (funts) Class
#'
#' The `funts` class is designed to encapsulate functional time series objects, including both univariate (FTS) and multivariate (MFTS) forms.
#' It provides a versatile framework for creating and manipulating `funts` objects, accommodating various basis systems and dimensions.
#'
#' @param X A matrix, three-dimensional array, or a list of matrix or array objects.
#'          When `method="data"`, it represents the observed curve values at discrete
#'          sampling points or argument values. When `method="coefs"`, `X` specifies the coefficients
#'          corresponding to the basis system defined in 'basisobj'. If `X` is a list, it defines a multivariate FTS,
#'          with each element being a matrix or three-dimensional array object. In matrix objects, rows correspond to argument values,
#'          and columns correspond to the length of the FTS. In three-dimensional array objects, the first and second dimensions correspond
#'          to argument values, and the third dimension to the length of the FTS.
#'
#' @param basisobj An object of class `basisfd`, a matrix of empirical basis, or a list of `basisfd` or empirical basis objects.
#' For empirical basis, rows correspond to basis functions, and columns correspond to grid points.
#'
#' @param argval A vector list of length `p` which is a set of argument values
#' corresponding to the observations in X. Each entry in this list should either
#' be a numeric value or a list of numeric elements, depending on the dimension of the domain
#' the variable is observed over. It can even vary from one variable to another,
#' If it be NULL, the default value for argval are the integers 1 to n, where n
#' is the size of the first dimension in argument X.?
#'
#' @param method Determines the type of the `X` matrix: "coefs" or "data."
#'
#' @param start The time of the first observation. It can be a single positive integer or an object of classes `Date`, `POSIXct`, or `POSIXt`,
#' specifying a natural time unit.
#'
#' @param end The time of the last observation, specified in the same way as `start`.
#' @param tname a string specifies the time index name
#' @param vnames a vector of strings specifies the variable names
#' @param dnames list of vector of strings specifies the variable domain names
#'
#'
#' @return An instance of the `funts` class containing functional time series data.
#'
#' @seealso \code{\link{fssa}}
#'
#' @examples
#' # 1D FTS example: Callcenter dataset
#' require(fda)
#' Call_data <- loadCallcenterData()
#' D <- matrix(sqrt(Call_data$calls), nrow = 240)
#' bs1 <- create.bspline.basis(c(0, 23), 22)
#' Y <- funts(X = D, basisobj = bs1, start = as.Date("1999-1-1"),
#'            vnames = "Sqrt of Call Numbers",
#'            dnames = "Time (6 minutes aggregated)",
#'            tname = "Date" )
#'
#' plot(Y, lwd = 2, npts = 200, col = "deepskyblue4",
#'   main = "Call Center Data"
#' )
#'
#' # _______2D Multivariate Example: Montana dataset
#' # Temperature curves and smoothed images of vegetation
#' montana_data <- loadMontanaData()
#' Temp <- montana_data$Temp
#' NDVI <- montana_data$NDVI
#' Montana_Data <- list(Temp / sd(Temp), NDVI)
#' bs1 <- create.bspline.basis(c(0, 23), 11)
#' bs2 <- create.bspline.basis(c(1, 33), 13)
#' bs2d <- list(bs2, bs2)
#' bsmv <- list(bs1, bs2d)
#' Y <- funts(X = Montana_Data, basisobj = bsmv,
#'            start = as.Date("2008-01-01"),
#'            end = as.Date("2013-09-30"),
#'            vnames = c("Normalized Temperature (\u00B0C)" , "NDVI"),
#'            dnames = list("Time", c("Latitude", "Longitude")),
#'            tname = "Date"
#' )
#'
#' plot(Y,  main = "Montana dataset")
#'
#' @return An instance of the `funts` class containing functional time series data.
#'
#' @export
funts <- function(X, basisobj, argval = NULL, method = "data", start = 1, end = NULL,
                  vnames = NULL, dnames = NULL, tname = NULL) { # Constructor function for the funts class

  # Check if X is a matrix, and if so, convert it to a list
  if (is.array(X)) {
    X <- list(X)
    basisobj <- list(basisobj)
    if (!is.null(argval)) argval <- list(argval)
  }
  if (is.basis(basisobj) | is.array(basisobj)) basisobj <- list(basisobj)
  if (is.numeric(argval)) argval <- list(argval)
  if (!is.list(X)) stop("The data input `X` must be a `matrix`,`array` or a list.")
  if (!is.list(basisobj)) stop("The basis input `basisobj` must be a `basisfd` object or a list of `basisfd` objects or empirical basis matrix.")
  # Check if all elements of X are matrices
  if (!all(sapply(X, is.array))) stop("All elements of the list X must be numeric `matrix` objects.")
  p <- length(X)
  dimSupp <- supp <- arg <- B_mat <- coefs <- list()
  n_def <- 100
  for (j in 1L:p) {
    # determine the dimension support of the variable j
    if (method == "data") {
      dimSupp[[j]] <- length(dim(X[[j]])) - 1
    } else if (is.list(argval[[j]])) {
      dimSupp[[j]] <- length(argval[[j]])
    } else if (!is.basis(basisobj[[j]]) && is.list(basisobj[[j]])) {
      dimSupp[[j]] <- length(basisobj[[j]])
    } else {
      dimSupp[[j]] <- 1
    }
    # Check vnames, dnames
    if (!is.null(vnames) & is.character(vnames) & length(vnames) != p) stop("The length of `vnames` must be equal to number of variables.")
    if (!is.null(dnames) & is.character(dnames) & length(dnames) == 1 & p == 1) dnames <- list(dnames)
    if (!is.null(dnames) & is.list(dnames) & length(dnames) != p) stop("The length of `dnames` must be equal to number of variables.")
    # Generating basis matrices=========================================
    # Generating a fd basis for variables whose domain is one-dimensional using a supplied list.
    if (dimSupp[[j]] == 1) { # 1d
      # setting up the grids:
      if (is.null(argval)) { # 1d with NULL grids
        if (is.basis(basisobj[[j]])) { # fd basis
          minval <- basisobj[[j]]$rangeval[1]
          maxval <- basisobj[[j]]$rangeval[2]
        } else { # Emp basis
          if (!is.null(attr(basisobj[[j]], "grids"))) {
            rangeval <- range(attr(basisobj[[j]], "grids"))
          } else if (!is.null(attr(basisobj[[j]], "rangeval"))) {
            rangeval <- attr(basisobj[[j]], "rangeval")
          } else {
            rangeval <- c(0, 1)
          }
          minval <- rangeval[1]
          maxval <- rangeval[2]
        }
        if (method == "data") { # method == "data" and NULL grids:
          arg[[j]] <- seq(from = minval, to = maxval, length.out = nrow(X[[j]]))
        } else { # method == "coefs" and NULL grids:
          arg[[j]] <- seq(from = minval, to = maxval, length.out = n_def)
        }
      } else { # 1d Not NULL grids:
        arg[[j]] <- argval[[j]]
      }
      if (is.basis(basisobj[[j]])) { # fd basis
        B_mat[[j]] <- eval.basis(evalarg = arg[[j]], basisobj = basisobj[[j]])
        dimnames(B_mat[[j]]) <- NULL
      } else {
        B_mat[[j]] <- eval.empb(evalarg = arg[[j]], basisobj = basisobj[[j]])[, ]
      }
    } else if (dimSupp[[j]] > 1) { # 2d
      # setting up u and v for grids:
      if (is.null(argval)) { # 2d with NULL grids
        if (all(sapply(basisobj[[j]], is.basis))) { # fd basis
          minval1 <- basisobj[[j]][[1]]$rangeval[1]
          maxval1 <- basisobj[[j]][[1]]$rangeval[2]
          minval2 <- basisobj[[j]][[2]]$rangeval[1]
          maxval2 <- basisobj[[j]][[2]]$rangeval[2]
        } else if (is.list(basisobj[[j]])) {
          if (!is.null(attr(basisobj[[j]][[1]], "grids"))) {
            rangeval1 <- range(attr(basisobj[[j]][[1]], "grids"))
          } else if (!is.null(attr(basisobj[[j]][[1]], "rangeval"))) {
            rangeval1 <- attr(basisobj[[j]][[1]], "rangeval")
          } else {
            rangeval1 <- c(0, 1)
          }
          if (!is.null(attr(basisobj[[j]][[2]], "grids"))) {
            rangeval2 <- range(attr(basisobj[[j]][[2]], "grids"))
          } else if (!is.null(attr(basisobj[[j]][[2]], "rangeval"))) {
            rangeval2 <- attr(basisobj[[j]][[2]], "rangeval")
          } else {
            rangeval2 <- c(0, 1)
          }
          minval1 <- rangeval1[1]
          maxval1 <- rangeval1[2]
          minval2 <- rangeval2[1]
          maxval2 <- rangeval2[2]
        } else {
          minval1 <- minval2 <- 0
          maxval1 <- maxval2 <- 1
        }
        if (method == "data") { # method == "data" and NULL grids:
          u <- seq(from = minval1, to = maxval1, length.out = dim(X[[j]])[1])
          v <- seq(from = minval2, to = maxval2, length.out = dim(X[[j]])[2])
        } else { # method == "coefs" and NULL grids:
          u <- seq(from = minval1, to = maxval1, length.out = n_def / 4)
          v <- seq(from = minval2, to = maxval2, length.out = n_def / 4)
        }
      } else { # 2d Not NULL grids:
        u <- argval[[j]][[1]]
        v <- argval[[j]][[2]]
      }
      if (all(sapply(basisobj[[j]], is.basis))) { # fd basis
        b_1 <- eval.basis(evalarg = u, basisobj = basisobj[[j]][[1]])
        b_2 <- eval.basis(evalarg = v, basisobj = basisobj[[j]][[2]])
      } else { # Empirical basis (list)
        if (is.list(basisobj[[j]])) {
          b_1 <- eval.empb(evalarg = u, basisobj = basisobj[[j]][[1]])
          b_2 <- eval.empb(evalarg = v, basisobj = basisobj[[j]][[2]])
        } else { # Empirical basis (Kronecker product or SVD)
          b_1 <- basisobj[[j]][, ]
          b_2 <- 1
        }
      }
      # Create the grid using u and v
      arg[[j]] <- cbind(rep(u, each = length(v)), rep(v, length(u)))
      B_mat[[j]] <- kronecker(b_1, b_2)
    }
    if (method == "data") {
      ### We should through an error if nrow(X), nrow(B) do not match in the empirical case ###
      # # Reshape funts of Images from arrays to matrices
      if (dimSupp[[j]] > 1) {
        M_x <- length(u)
        M_y <- length(v)
        if (is.matrix(X[[j]])) X[[j]] <- array(X[[j]], dim = c(M_x, M_y, 1))
        X[[j]] <- matrix(aperm(X[[j]], c(2, 1, 3)), nrow = M_x * M_y)
      }
      # Estimate the coefficients of each funts variables.=========================================
      lambda_min_gram = check_basis(B_mat[[j]])
      if(lambda_min_gram < 10^-12){
        stop(paste0("The provided basis for variable ",as.character(j)," is not linearly independent."))
      }
      coefs[[j]] <- solve(t(B_mat[[j]]) %*% B_mat[[j]]) %*% t(B_mat[[j]]) %*% X[[j]]
    } else { # method == "coefs"
      coefs[[j]] <- X[[j]]
    }
  }
  # Creating the time object ========================================
  N <- tail(dim(X[[1]]), 1)
  if (is.null(end)) end <- start + N - 1
  time <- seq(from = start, to = end, length.out = N)

  # Create and return an instance of the funts class=========================================
  out <- list(N = N, dimSupp = dimSupp, time = time, coefs = coefs, basis = basisobj,
              B_mat = B_mat, argval = arg, vnames = vnames, dnames = dnames, tname = tname)
  class(out) <- "funts"
  return(out)
}
