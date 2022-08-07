#' Functional Time Series Class
#'
#' This function is used to create functional time series (\code{\link{fts}}) objects from lists of discretely sampled data, basis specifications, and grid elements which provide the domain that each variable is observed over. Each variable is assumed to be observed over a regular and equidistant grid. In addition, each variable in the fts is assumed to be observed over a one or two-dimensional domain.
#' @param X  A list of length p where p is a positive integer specifying the number of variables observed in the fts. Each entry in the list should be a matrix or an array.
#' @param B A list of length p. Each entry in the list should be either a matrix specifying the basis for each variable or each list entry should be a list specifying the number of basis elements and desired basis type to be used in the smoothing process.
#' @param grid A list of length p. Each entry in the list should either be a numeric or a list of numeric elements depending on the dimension of the domain the variable is observed over.
#' @param time A character vector where each entry specifies the time an observation is made.
#' @importFrom fda create.bspline.basis create.fourier.basis eval.basis
#' @importFrom methods is new
#' @seealso \code{\link{fssa}}
#' @note Refer to \code{\link{fssa}} for examples on how to run this function.
#' @export
fts <- function(X, B, grid, time = NULL) {
  fts <- new(Class = "fts", C = X, B = B, grid = grid, basis_type = B)
  p <- length(X)
  # Loop over each variable in the fts
  for (j in 1L:p) {
    if(isFALSE(is.list(B[[j]]))) fts@basis_type[[j]] = "User-defined basis.";
    # Deciding whether or not the variable corresponds to an fts observed over a two-dimensional domain or a one-dimensional domain.
    if (is.list(grid[[j]]) && length(grid[[j]]) == 2) {
      if (length(grid[[j]][[1]]) == 2 && length(grid[[j]][[2]]) == 2) {
        # If the grid is specified with min and max values in both horizontal and vertical axes.
        # Example: list(c(1,2),c(3,4))

        u <- seq(grid[[j]][[1]][1], grid[[j]][[1]][2], length.out = ncol(X[[j]][,,1]))
        v <- seq(grid[[j]][[2]][1], grid[[j]][[2]][2], length.out = nrow(X[[j]][,,1]))
      } else if (length(grid[[j]][[1]]) == 2 && length(grid[[j]][[2]]) > 2) {
        # If the grid is specified with min and max values in the horizontal axis and with the actual grid in the vertical axis.
        # Example: list(c(1,2),v)

        u <- seq(grid[[j]][[1]][1], grid[[j]][[1]][2], length.out = ncol(X[[j]][,,1]))
        v <- grid[[j]][[2]]
      } else if (length(grid[[j]][[1]]) > 2 && length(grid[[j]][[2]]) == 2) {
        # This case is similar as seen in the previous else if statement
        # Example: list(u,c(1,2))

        u <- grid[[j]][[1]]
        v <- seq(grid[[j]][[2]][1], grid[[j]][[2]][2], length.out = nrow(X[[j]][,,1]))
      } else {

        # If the grid is specified with the actual grid in the horizontal and vertical axes.
        # Example: list(u,v)
        u <- seq(min(grid[[j]][[1]]), max(grid[[j]][[1]]), length.out = length(grid[[j]][[1]]))
        v <- seq(min(grid[[j]][[2]]), max(grid[[j]][[2]]), length.out = length(grid[[j]][[2]]))
      }
      M_x <- length(u)
      M_y <- length(v)

      # Create the grid using u and v
      fts@grid[[j]] <- cbind(rep(u, each = length(v)), rep(v, length(u)))

      # Reshape FTS of Images from matricies to vectors
      if (is.matrix(X[[j]])) X[[j]] <- array(X[[j]], dim = c(M_x, M_y, 1))
      X[[j]] <- matrix(aperm(X[[j]], c(2, 1, 3)), nrow = M_x * M_y)

      # If the variable is observed over a one-dimensional domain
    } else if (is.list(grid[[j]]) && length(grid[[j]]) == 1 && length(grid[[j]][[1]]) > 2) {
      # If the grid is supplied in a list, the variable corresponds with a one-dimensional domain, and the list element is the grid.
      # Example: list(u)
      M_x <- length(grid[[j]][[1]])
      u <- seq(min(grid[[j]][[1]]), max(grid[[j]][[1]]), length.out = M_x)
      fts@grid[[j]] <- matrix(u, ncol = 1)
      X[[j]] <- as.matrix(X[[j]])
    } else if (is.list(grid[[j]]) && length(grid[[j]]) == 1 && length(grid[[j]][[1]]) == 2) {

      # If the grid is supplied in a list, the variable corresponds with a one-dimensional domain, and the list element is the min and max values of the grid.
      # Example: list(c(1,2))
      u <- seq(grid[[j]][[1]][1], grid[[j]][[1]][2], length.out = nrow(X[[j]]))
      fts@grid[[j]] <- matrix(u, ncol = 1)
      X[[j]] <- as.matrix(X[[j]])
    } else if (is.numeric(grid[[j]]) && length(grid[[j]]) == 2) {

      # If the grid specified with the min and max values of the grid.
      # Example: c(1,2)
      u <- seq(grid[[j]][1], grid[[j]][2], length.out = nrow(X[[j]]))
      fts@grid[[j]] <- matrix(u, ncol = 1)
      X[[j]] <- as.matrix(X[[j]])
    } else {
      # If the grid specified with the grid itself.
      # Example: u
      fts@grid[[j]] <- matrix(grid[[j]], ncol = 1)
      X[[j]] <- as.matrix(X[[j]])
    }



    # Generating basis matrices
    if (is.list(B[[j]]) && B[[j]][[2]] == "bspline" && ncol(fts@grid[[j]]) == 1) {
      # Generating a bspline basis for variables whose domain is one-dimensional using a supplied list.
      # Example: list(10,"bspline")
      fts@B[[j]] <- fda::eval.basis(evalarg = seq(min(fts@grid[[j]]), max(fts@grid[[j]]), length.out = nrow(fts@grid[[j]])), basisobj = fda::create.bspline.basis(rangeval = c(min(fts@grid[[j]]), max(fts@grid[[j]])), nbasis = B[[j]][[1]]))
      B[[j]] <- fts@B[[j]]
    } else if (is.list(B[[j]]) && B[[j]][[2]] == "fourier" && ncol(fts@grid[[j]]) == 1) {
      # Generating a fourier basis for variables whose domain is one-dimensional using a supplied list.
      # Example: list(10,"fourier")
      fts@B[[j]] <- fda::eval.basis(evalarg = seq(min(fts@grid[[j]]), max(fts@grid[[j]]), length.out = nrow(fts@grid[[j]])), basisobj = fda::create.fourier.basis(rangeval = c(min(fts@grid[[j]]), max(fts@grid[[j]])), nbasis = B[[j]][[1]]))
      B[[j]] <- fts@B[[j]]
    }

    if (is.list(B[[j]]) && B[[j]][[3]] == "bspline" && B[[j]][[4]] == "bspline" && ncol(fts@grid[[j]]) == 2) {
      # Generating a bspline basis for variables whose domain is two-dimensional using a supplied list.
      # Example: list(list(10,"bspline"),list(12,"bspline"))
      b_1 <- fda::eval.basis(evalarg = seq(min(fts@grid[[j]][, 1]), max(fts@grid[[j]][, 1]), length.out = M_x), basisobj = fda::create.bspline.basis(rangeval = c(min(fts@grid[[j]][, 1]), max(fts@grid[[j]][, 1])), nbasis = B[[j]][[1]]))
      b_2 <- fda::eval.basis(evalarg = seq(min(fts@grid[[j]][, 2]), max(fts@grid[[j]][, 2]), length.out = M_y), basisobj = fda::create.bspline.basis(rangeval = c(min(fts@grid[[j]][, 2]), max(fts@grid[[j]][, 2])), nbasis = B[[j]][[2]]))
      fts@B[[j]] <- kronecker(b_1, b_2)
      B[[j]] <- fts@B[[j]]
    } else if (is.list(B[[j]]) && B[[j]][[3]] == "bspline" && B[[j]][[4]] == "fourier" && ncol(fts@grid[[j]]) == 2) {
      # Generating a bspline/fourier basis for variables whose domain is two-dimensional using a supplied list.
      # Example: list(list(10,"bspline"),list(12,"fourier"))
      b_1 <- fda::eval.basis(evalarg = seq(min(fts@grid[[j]][, 1]), max(fts@grid[[j]][, 1]), length.out = M_x), basisobj = fda::create.bspline.basis(rangeval = c(min(fts@grid[[j]][, 1]), max(fts@grid[[j]][, 1])), nbasis = B[[j]][[1]]))
      b_2 <- fda::eval.basis(evalarg = seq(min(fts@grid[[j]][, 2]), max(fts@grid[[j]][, 2]), length.out = M_y), basisobj = fda::create.fourier.basis(rangeval = c(min(fts@grid[[j]][, 2]), max(fts@grid[[j]][, 2])), nbasis = B[[j]][[2]]))
      fts@B[[j]] <- kronecker(b_1, b_2)
      B[[j]] <- fts@B[[j]]
    } else if (is.list(B[[j]]) && B[[j]][[3]] == "fourier" && B[[j]][[4]] == "bspline" && ncol(fts@grid[[j]]) == 2) {
      # This case is similar to the previous else if statement
      # Example: list(list(10,"fourier"),list(12,"bspline"))
      b_1 <- fda::eval.basis(evalarg = seq(min(fts@grid[[j]][, 1]), max(fts@grid[[j]][, 1]), length.out = M_x), basisobj = fda::create.fourier.basis(rangeval = c(min(fts@grid[[j]][, 1]), max(fts@grid[[j]][, 1])), nbasis = B[[j]][[1]]))
      b_2 <- fda::eval.basis(evalarg = seq(min(fts@grid[[j]][, 2]), max(fts@grid[[j]][, 2]), length.out = M_y), basisobj = fda::create.bspline.basis(rangeval = c(min(fts@grid[[j]][, 2]), max(fts@grid[[j]][, 2])), nbasis = B[[j]][[2]]))
      fts@B[[j]] <- kronecker(b_1, b_2)
      B[[j]] <- fts@B[[j]]
    } else if (is.list(B[[j]]) && B[[j]][[3]] == "fourier" && B[[j]][[4]] == "fourier" && ncol(fts@grid[[j]]) == 2) {
      # Generating a fourier basis for variables whose domain is two-dimensional using a supplied list.
      # Example: list(list(10,"fourier"),list(12,"fourier"))
      b_1 <- fda::eval.basis(evalarg = seq(min(fts@grid[[j]][, 1]), max(fts@grid[[j]][, 1]), length.out = M_x), basisobj = fda::create.fourier.basis(rangeval = c(min(fts@grid[[j]][, 1]), max(fts@grid[[j]][, 1])), nbasis = B[[j]][[1]]))
      b_2 <- fda::eval.basis(evalarg = seq(min(fts@grid[[j]][, 2]), max(fts@grid[[j]][, 2]), length.out = M_y), basisobj = fda::create.fourier.basis(rangeval = c(min(fts@grid[[j]][, 2]), max(fts@grid[[j]][, 2])), nbasis = B[[j]][[2]]))
      fts@B[[j]] <- kronecker(b_1, b_2)
      B[[j]] <- fts@B[[j]]
    }
    # Estimate the coefficients of each fts variable using basis elements and observed data.

    fts@C[[j]] <- solve(t(B[[j]]) %*% B[[j]]) %*% t(B[[j]]) %*% X[[j]]
    if(is.null(time) && is.null(colnames(fts@C[[j]]))){
      time = as.character(1:ncol(fts@C[[j]]))
    }else if(is.null(time) && is.null(colnames(fts@C[[j]]))==FALSE){
      time = colnames(fts@C[[j]])
    }
    colnames(fts@C[[j]]) = as.character(time)

    if (is.null(rownames(fts@grid[[j]])) == TRUE) {
      rownames(fts@grid[[j]]) <- as.character(1:nrow(fts@grid[[j]]))
    }
  }


  return(fts)
}
# class definitions, methods, etc. go outside of constructor.

check_fts <- function(object) {
  errors <- character()
  msg <- NULL

  # check if fts lengths are equal
  if (length(object@C) > 1) {
    for (j_1 in 1:length(object@C)) {
      for (j_2 in j_1:length(object@C)) {
        if (length(dim(object@C[[j_1]])) == 2 && length(dim(object@C[[j_2]])) == 2) {
          if (ncol(object@C[[j_1]]) != ncol(object@C[[j_2]])) {
            msg <- paste("The length of variable", as.character(j_1), "and the length of variable", as.character(j_2), "are not equal.")

            errors <- c(errors, msg)

            msg <- NULL
          }
        } else if (length(dim(object@C[[j_1]])) == 2 && length(dim(object@C[[j_2]])) == 3) {
          if (ncol(object@C[[j_1]]) != dim(object@C[[j_2]])[3]) {
            msg <- paste("The length of variable", as.character(j_1), "and the length of variable", as.character(j_2), "are not equal.")

            errors <- c(errors, msg)

            msg <- NULL
          }
        } else if (length(dim(object@C[[j_1]])) == 3 && length(dim(object@C[[j_2]])) == 2) {
          if (dim(object@C[[j_1]])[3] != ncol(object@C[[j_2]])) {
            msg <- paste("The length of variable", as.character(j_1), "and the length of variable", as.character(j_2), "are not equal.")

            errors <- c(errors, msg)

            msg <- NULL
          }
        } else if (length(dim(object@C[[j_1]])) == 3 && length(dim(object@C[[j_2]])) == 3) {
          if (dim(object@C[[j_1]])[3] != dim(object@C[[j_2]])[3]) {
            msg <- paste("The length of variable", as.character(j_1), "and the length of variable", as.character(j_2), "are not equal.")

            errors <- c(errors, msg)

            msg <- NULL
          }
        }
      }
    }
  }

  # check if provided lists have same length
  if (length(object@C) != length(object@B)) {
    msg <- paste("The length of X is", as.character(length(object@C)), "and the length of B is", as.character(length(object@B)), "whereas the length of these lists should be equal.")

    errors <- c(errors, msg)

    msg <- NULL
  }
  if (length(object@C) != length(object@grid)) {
    msg <- paste("The length of X is", as.character(length(object@C)), "and the length of grid is", as.character(length(object@grid)), "whereas the length of these lists should be equal.")

    errors <- c(errors, msg)

    msg <- NULL
  }

  if (length(object@B) != length(object@grid)) {
    msg <- paste("The length of B is", as.character(length(object@B)), "and the length of grid is", as.character(length(object@grid)), "whereas the length of these lists should be equal.")

    errors <- c(errors, msg)

    msg <- NULL
  }

  # check if the classes and values for supplied inputs are correct
  for (i in 1L:length(object@C)) {
    if (isFALSE(is.matrix(object@C[[i]])) && isFALSE(is.numeric(object@C[[i]])) && isFALSE(is.integer(object@C[[i]])) && isFALSE(is.array(object@C[[i]]))) {
      msg <- paste("The class of entry", as.character(i), "of X is not matrix, numeric, or integer as is expected.")

      errors <- c(errors, msg)

      msg <- NULL
    }
  }

  for (i in 1L:length(object@B)) {
    if (isFALSE(is.matrix(object@B[[i]])) && isFALSE(is.numeric(object@B[[i]])) && isFALSE(is.integer(object@B[[i]])) && isFALSE(is.array(object@B[[i]])) && isFALSE(is.list(object@B[[i]]))) {
      msg <- paste("The class of entry", as.character(i), "of B is not matrix, numeric, integer, or list as is expected.")

      errors <- c(errors, msg)

      msg <- NULL
    }

    if (is.list(object@grid[[i]]) && length(object@grid[[i]]) == 2) {
      if (is.list(object@B[[i]]) && length(object@B[[i]]) != 4) {
        msg <- paste("The length of entry", as.character(i), "of B is not four. The entry should be a basis matrix or a list of the form list(int_1,int_2,string_1,string_2) for variables taken over a two-dimensional domain where int_1,int_2>0 and string_1,string_2 = \"bspline\" or string_1,string_2 = \"fourier\".")

        errors <- c(errors, msg)

        msg <- NULL
      }

      if (is.list(object@B[[i]]) && length(object@B[[i]]) == 4 && isFALSE(is.integer(object@B[[i]][[1]])) && isFALSE(is.numeric(object@B[[i]][[1]]))) {
        msg <- paste("The first element in entry", as.character(i), "of B should be a positive integer.")

        errors <- c(errors, msg)

        msg <- NULL
      }

      if ((is.list(object@B[[i]]) && length(object@B[[i]]) == 4 && is.integer(object@B[[i]][[1]]) && object@B[[i]][[1]] <= 0) || (is.list(object@B[[i]]) && length(object@B[[i]]) == 4 && is.numeric(object@B[[i]][[1]]) && object@B[[i]][[1]] <= 0)) {
        msg <- paste("The first element in entry", as.character(i), "of B should be a positive integer.")

        errors <- c(errors, msg)

        msg <- NULL
      }


      if (is.list(object@B[[i]]) && length(object@B[[i]]) == 4 && isFALSE(is.integer(object@B[[i]][[2]])) && isFALSE(is.numeric(object@B[[i]][[2]]))) {
        msg <- paste("The second element in entry", as.character(i), "of B should be a positive integer.")

        errors <- c(errors, msg)

        msg <- NULL
      }

      if ((is.list(object@B[[i]]) && length(object@B[[i]]) == 4 && is.integer(object@B[[i]][[2]]) && object@B[[i]][[2]] <= 0) || (is.list(object@B[[i]]) && length(object@B[[i]]) == 4 && is.numeric(object@B[[i]][[2]]) && object@B[[i]][[2]] <= 0)) {
        msg <- paste("The second element in entry", as.character(i), "of B should be a positive integer.")

        errors <- c(errors, msg)

        msg <- NULL
      }

      if (is.list(object@B[[i]]) && length(object@B[[i]]) == 4 && isFALSE(is.character(object@B[[i]][[3]]))) {
        msg <- paste("The third element in entry", as.character(i), "of B should be a string that is either \"bspline\" or \"fourier\".")

        errors <- c(errors, msg)

        msg <- NULL
      }

      if (is.list(object@B[[i]]) && length(object@B[[i]]) == 4 && is.character(object@B[[i]][[3]]) && object@B[[i]][[3]] != "bspline" && object@B[[i]][[3]] != "fourier") {
        msg <- paste("The third element in entry", as.character(i), "of B should be a string that is either \"bspline\" or \"fourier\".")

        errors <- c(errors, msg)

        msg <- NULL
      }

      if (is.list(object@B[[i]]) && length(object@B[[i]]) == 4 && isFALSE(is.character(object@B[[i]][[4]]))) {
        msg <- paste("The fourth element in entry", as.character(i), "of B should be a string that is either \"bspline\" or \"fourier\".")

        errors <- c(errors, msg)

        msg <- NULL
      }

      if (is.list(object@B[[i]]) && length(object@B[[i]]) == 4 && is.character(object@B[[i]][[4]]) && object@B[[i]][[4]] != "bspline" && object@B[[i]][[4]] != "fourier") {
        msg <- paste("The fourth element in entry", as.character(i), "of B should be a string that is either \"bspline\" or \"fourier\".")

        errors <- c(errors, msg)

        msg <- NULL
      }
    } else {
      if (is.list(object@B[[i]])[[1]] && length(object@B[[i]]) != 2) {
        msg <- paste("The length of entry", as.character(i), "of B is not two. The entry should be a basis matrix or a list of the form list(int,string) for variables taken over a one-dimensional domain where int>0 and string = \"bspline\" or string = \"fourier\".")

        errors <- c(errors, msg)

        msg <- NULL
      }

      if (is.list(object@B[[i]]) && length(object@B[[i]]) == 2 && isFALSE(is.integer(object@B[[i]][[1]])) && isFALSE(is.numeric(object@B[[i]][[1]]))) {
        msg <- paste("The first element in entry", as.character(i), "of B should be a positive integer.")

        errors <- c(errors, msg)

        msg <- NULL
      }

      if ((is.list(object@B[[i]]) && length(object@B[[i]]) == 2 && is.integer(object@B[[i]][[1]]) && object@B[[i]][[1]] <= 0) || (is.list(object@B[[i]]) && length(object@B[[i]]) == 2 && is.numeric(object@B[[i]][[1]]) && object@B[[i]][[1]] <= 0)) {
        msg <- paste("The first element in entry", as.character(i), "of B should be a positive integer.")

        errors <- c(errors, msg)

        msg <- NULL
      }

      if (is.list(object@B[[i]]) && length(object@B[[i]]) == 2 && isFALSE(is.character(object@B[[i]][[2]]))) {
        msg <- paste("The second element in entry", as.character(i), "of B should be a string that is either \"bspline\" or \"fourier\".")

        errors <- c(errors, msg)

        msg <- NULL
      }

      if (is.list(object@B[[i]]) && length(object@B[[i]]) == 2 && is.character(object@B[[i]][[2]]) && object@B[[i]][[2]] != "bspline" && object@B[[i]][[2]] != "fourier") {
        msg <- paste("The second element in entry", as.character(i), "of B should be a string that is either \"bspline\" or \"fourier\".")

        errors <- c(errors, msg)

        msg <- NULL
      }
    }
  }
  # Check classes and supplied inputs for the grids
  for (i in 1L:length(object@grid)) {
    if (is.list(object@grid[[i]])) {
      for (k in 1:length(object@grid[[i]])) {
        if (isFALSE(is.matrix(object@grid[[i]][[k]])) && isFALSE(is.integer(object@grid[[i]][[k]])) && isFALSE(is.numeric(object@grid[[i]][[k]]))) {
          msg <- paste("Element ", as.character(k), " of entry ", as.character(i), " of grid is not a matrix, numeric, or integer as is expected.", sep = "")
          errors <- c(errors, msg)

          msg <- NULL
        }

        if (is.numeric(object@grid[[i]][[k]]) && object@grid[[i]][[k]][1] >= object@grid[[i]][[k]][2]) {
          msg <- paste("The first element of interval ", as.character(k), " of entry ", as.character(i), " of grid should be less than the second element of that same interval.", sep = "")

          errors <- c(errors, msg)

          msg <- NULL
        }
      }
    } else {
      if (isFALSE(is.matrix(object@grid[[i]])) && isFALSE(is.numeric(object@grid[[i]])) && isFALSE(is.integer(object@grid[[i]])[[1]])) {
        msg <- paste("The class of entry ", as.character(i), " of grid is not matrix, numeric, integer, or list as is expected.", sep = "")

        errors <- c(errors, msg)

        msg <- NULL
      }


      if (is.numeric(object@grid[[i]]) && object@grid[[i]][1] >= object@grid[[i]][2]) {
        msg <- paste("The first element of entry ", as.character(i), " of grid should be less than the second element.", sep = "")

        errors <- c(errors, msg)

        msg <- NULL
      }
    }
  }

  if (length(errors) != 0) {
    return(errors)
  }
}

# Set the class definition
setClass("fts", representation(C = "list", B = "list", grid = "list", basis_type = "list"), validity = check_fts)
# Show method for the fts object
setMethod("show", "fts", function(object) {
  N <- ncol(object@C[[1]])
  p <- length(object@C)

  if (length(object@C) == 1) {
    cat("Univariate ", is(object), " of length ", N, ".", "\n", sep = "")
  } else {
    cat("Multivariate ", is(object), " of length ", N, " that is comprised of ", p, " variables.", "\n", sep = "")
  }

  for (j in 1:p) {
    if (ncol(object@grid[[j]]) == 1) {

      object@grid[[j]] <- matrix(data = object@grid[[j]], ncol = 1)
      number_points <- nrow(object@grid[[j]])
      min_x <- object@grid[[j]][1, 1]
      max_x <- object@grid[[j]][number_points, 1]

      cat("\n", "The domain of the observations for variable ", as.character(j), " is one-dimensional.", "\n",
        "Observations for variable ", as.character(j), " are sampled on a regular grid at ", as.character(number_points), " points.", "\n",
        "The x axis for this variable varies along the domain [", as.character(min_x), ",", as.character(max_x), "].", "\n",
        sep = ""
      )
    } else {

      number_points <- nrow(object@grid[[j]])
      min_x <- object@grid[[j]][1, 1]
      min_y <- object@grid[[j]][1, 2]
      max_x <- object@grid[[j]][number_points, 1]
      max_y <- object@grid[[j]][number_points, 2]

      cat("\n", "The domain of the observations for variable ", as.character(j), " is two-dimensional.", "\n",
        "Observations for variable ", as.character(j), " are sampled on a regular grid at ", as.character(number_points), " points.", "\n",
        "The x axis for this variable varies along the domain [", as.character(min_x), ",", as.character(max_x), "].", "\n",
        "The y axis for this variable varies along the domain [", as.character(min_y), ",", as.character(max_y), "].", "\n",
        sep = ""
      )
    }
  }
})
