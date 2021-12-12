#' Addition of Functional Time Series
#'
#' A method for functional time series (\code{\link{fts}}) addition and fts-scalar addition. Note that if the fts is multivariate
#' then a vector of numerics may be provided allowing for addition of different scalars to different variables. For example, multivariate fts-numeric addition
#' follows the form of \code{Y+c(1,2)} if \code{Y} is a bivariate fts.
#' @return An object of class \code{\link{fts}}.
#'
#' @param Y1 An object of class \code{\link{fts}} or numeric.
#' @param Y2 An object of class \code{\link{fts}} or numeric.
#' @examples
#' \dontrun{
#' require(Rfssa)
#' load_github_data("https://github.com/haghbinh/Rfssa/blob/main/data/Callcenter.RData")
#' D <- matrix(sqrt(Callcenter$calls), nrow = 240)
#' u <- seq(0, 1, length.out = 240) # Define domain of functional data
#' d <- 22 # number of basis elements
#' Y <- fts(list(D), list(list(d, "bspline")), list(u))
#' plot(Y)
#' Yplus <- Y + Y # add the functional time series to itself
#' plot(Yplus)
#' Yplus2 <- Y + 2 # add 2 to every term in the functional time series
#' plot(Yplus2)
#' }
#' @seealso \code{\link{fts}}
#' @useDynLib Rfssa
#' @export
"+.fts" <- function(Y1, Y2) {
  if (class(Y1) == "fts" && class(Y2) == "fts") {
    if (length(Y1@C) != length(Y2@C)) {
      stop("Functional time series are of different length.")
    }


    p <- length(Y1@C)
    grid <- Y1@grid
    basis <- Y1@B
    eval_fts <- list()
    new_grid <- list()
    for (j in 1:p) {
      if (max(max(Y1@B[[j]] - Y2@B[[j]])) >= 10^-14) {
        stop(paste("The basis matrix corresponding to variable ", as.character(j), " should be the same between fts objects.", sep = ""))
      }

      eval_fts[[j]] <- basis[[j]] %*% (Y1@C[[j]] + Y2@C[[j]])
      N <- ncol(eval_fts[[j]])

      if (ncol(Y1@grid[[j]]) == 2) {
        x <- unique(Y1@grid[[j]][, 1])
        y <- unique(Y1@grid[[j]][, 2])
        new_fts_two_d <- array(data = NA, dim = c(length(x), length(y), N))
        for (n in 1:N) {
          count <- 1
          for (i_1 in 1:length(x)) {
            for (i_2 in 1:length(y)) {
              new_fts_two_d[i_1, i_2, n] <- eval_fts[[j]][count, n]
              count <- count + 1
            }
          }
        }

        eval_fts[[j]] <- new_fts_two_d
        new_grid[[j]] <- list(x, y)
      } else {
        new_grid[[j]] <- Y1@grid[[j]]
      }
    }
  } else if (class(Y1) == "fts" && class(Y2) == "numeric" || class(Y1) == "fts" && class(Y2) == "matrix") {
    Y2 <- as.numeric(Y2)
    p <- length(Y1@C)
    if (length(Y2) == 1) Y2 <- rep(Y2, p)
    grid <- Y1@grid
    basis <- Y1@B
    eval_fts <- list()
    new_grid <- list()
    for (j in 1:p) {
      eval_fts[[j]] <- basis[[j]] %*% (Y1@C[[j]]) + Y2[j]
      N <- ncol(eval_fts[[j]])

      if (ncol(Y1@grid[[j]]) == 2) {
        x <- unique(Y1@grid[[j]][, 1])
        y <- unique(Y1@grid[[j]][, 2])
        new_fts_two_d <- array(data = NA, dim = c(length(x), length(y), N))
        for (n in 1:N) {
          count <- 1
          for (i_1 in 1:length(x)) {
            for (i_2 in 1:length(y)) {
              new_fts_two_d[i_1, i_2, n] <- eval_fts[[j]][count, n]
              count <- count + 1
            }
          }
        }

        eval_fts[[j]] <- new_fts_two_d
        new_grid[[j]] <- list(x, y)
      } else {
        new_grid[[j]] <- Y1@grid[[j]]
      }
    }
  } else if (class(Y1) == "numeric" && class(Y2) == "fts" || class(Y1) == "matrix" && class(Y2) == "fts") {
    Y1 <- as.numeric(Y1)
    p <- length(Y2@C)
    if (length(Y1) == 1) Y1 <- rep(Y1, p)
    grid <- Y2@grid
    basis <- Y2@B
    eval_fts <- list()
    new_grid <- list()
    for (j in 1:p) {
      eval_fts[[j]] <- basis[[j]] %*% (Y2@C[[j]]) + Y1[j]
      N <- ncol(eval_fts[[j]])
      if (ncol(Y2@grid[[j]]) == 2) {
        x <- unique(Y2@grid[[j]][, 1])
        y <- unique(Y2@grid[[j]][, 2])
        new_fts_two_d <- array(data = NA, dim = c(length(x), length(y), N))
        for (n in 1:N) {
          count <- 1
          for (i_1 in 1:length(x)) {
            for (i_2 in 1:length(y)) {
              new_fts_two_d[i_1, i_2, n] <- eval_fts[[j]][count, n]
              count <- count + 1
            }
          }
        }

        eval_fts[[j]] <- new_fts_two_d
        new_grid[[j]] <- list(x, y)
      } else {
        new_grid[[j]] <- Y2@grid[[j]]
      }
    }
  }

  out <- Rfssa::fts(eval_fts, basis, new_grid)
  return(out)
}

#' Subtraction of Functional Time Series
#'
#' A method for functional time series (\code{\link{fts}}) subtraction and fts-scalar subtraction. Note that if the fts is multivariate
#' then a vector of numerics may be provided allowing for subtraction of different scalars from different variables. For example, multivariate fts-numeric subtraction
#' follows the form of \code{Y-c(1,2)} if \code{Y} is a bivariate fts.
#' @return An object of class \code{\link{fts}}.
#'
#'
#' @param Y1 An object of class \code{\link{fts}} or numeric.
#' @param Y2 An object of class \code{\link{fts}} or numeric.
#' @examples
#' \dontrun{
#' require(Rfssa)
#' load_github_data("https://github.com/haghbinh/Rfssa/blob/main/data/Callcenter.RData")
#' D <- matrix(sqrt(Callcenter$calls), nrow = 240)
#' u <- seq(0, 1, length.out = 240) # Define domain of functional data
#' d <- 22 # number of basis elements
#' Y <- fts(list(D), list(list(d, "bspline")), list(u))
#' plot(Y)
#' Yminus <- Y[1:4] - Y[5:8] # subtract the functional time series to itself
#' plot(Yminus)
#' Yminus2 <- Y - 2 # subtract 2 to every term in the functional time series
#' plot(Yminus2)
#' }
#'
#' @seealso \code{\link{fts}}
#' @useDynLib Rfssa
#' @export
"-.fts" <- function(Y1, Y2) {
  if (class(Y1) == "fts" && class(Y2) == "fts") {
    if (length(Y1@C) != length(Y2@C)) {
      stop("Functional time series are of different length.")
    }

    p <- length(Y1@C)
    grid <- Y1@grid
    basis <- Y1@B
    eval_fts <- list()
    new_grid <- list()
    for (j in 1:p) {
      if (max(max(Y1@B[[j]] - Y2@B[[j]])) >= 10^-14) {
        stop(paste("The basis matrix corresponding to variable ", as.character(j), " should be the same between fts objects.", sep = ""))
      }

      eval_fts[[j]] <- basis[[j]] %*% (Y1@C[[j]] - Y2@C[[j]])
      N <- ncol(eval_fts[[j]])
      if (ncol(Y1@grid[[j]]) == 2) {
        x <- unique(Y1@grid[[j]][, 1])
        y <- unique(Y1@grid[[j]][, 2])
        new_fts_two_d <- array(data = NA, dim = c(length(x), length(y), N))
        for (n in 1:N) {
          count <- 1
          for (i_1 in 1:length(x)) {
            for (i_2 in 1:length(y)) {
              new_fts_two_d[i_1, i_2, n] <- eval_fts[[j]][count, n]
              count <- count + 1
            }
          }
        }

        eval_fts[[j]] <- new_fts_two_d
        new_grid[[j]] <- list(x, y)
      } else {
        new_grid[[j]] <- Y1@grid[[j]]
      }
    }
  } else if (class(Y1) == "fts" && class(Y2) == "numeric" || class(Y1) == "fts" && class(Y2) == "matrix") {
    Y2 <- as.numeric(Y2)
    p <- length(Y1@C)
    if (length(Y2) == 1) Y2 <- rep(Y2, p)
    grid <- Y1@grid
    basis <- Y1@B
    eval_fts <- list()
    new_grid <- list()
    for (j in 1:p) {
      eval_fts[[j]] <- basis[[j]] %*% (Y1@C[[j]]) - Y2[j]
      N <- ncol(eval_fts[[j]])

      if (ncol(Y1@grid[[j]]) == 2) {
        x <- unique(Y1@grid[[j]][, 1])
        y <- unique(Y1@grid[[j]][, 2])
        new_fts_two_d <- array(data = NA, dim = c(length(x), length(y), N))
        for (n in 1:N) {
          count <- 1
          for (i_1 in 1:length(x)) {
            for (i_2 in 1:length(y)) {
              new_fts_two_d[i_1, i_2, n] <- eval_fts[[j]][count, n]
              count <- count + 1
            }
          }
        }

        eval_fts[[j]] <- new_fts_two_d
        new_grid[[j]] <- list(x, y)
      } else {
        new_grid[[j]] <- Y1@grid[[j]]
      }
    }
  } else if (class(Y1) == "numeric" && class(Y2) == "fts" || class(Y1) == "matrix" && class(Y2) == "fts") {
    Y1 <- as.numeric(Y1)
    p <- length(Y2@C)
    if (length(Y1) == 1) Y1 <- rep(Y1, p)
    grid <- Y2@grid
    basis <- Y2@B
    eval_fts <- list()
    new_grid <- list()
    for (j in 1:p) {
      eval_fts[[j]] <- Y1[j] - basis[[j]] %*% (Y2@C[[j]])
      N <- ncol(eval_fts[[j]])

      if (ncol(Y2@grid[[j]]) == 2) {
        x <- unique(Y2@grid[[j]][, 1])
        y <- unique(Y2@grid[[j]][, 2])
        new_fts_two_d <- array(data = NA, dim = c(length(x), length(y), N))
        for (n in 1:N) {
          count <- 1
          for (i_1 in 1:length(x)) {
            for (i_2 in 1:length(y)) {
              new_fts_two_d[i_1, i_2, n] <- eval_fts[[j]][count, n]
              count <- count + 1
            }
          }
        }

        eval_fts[[j]] <- new_fts_two_d
        new_grid[[j]] <- list(x, y)
      } else {
        new_grid[[j]] <- Y2@grid[[j]]
      }
    }
  }

  out <- Rfssa::fts(eval_fts, basis, new_grid)
  return(out)
}


#' Multiplication of Functional Time Series
#'
#' A method for functional time series (\code{\link{fts}}) pointwise multiplication and fts-scalar multiplication. Note that if the fts is multivariate
#' then a vector of numerics may be provided allowing for multiplication of different variables by different scalars. For example, multivariate fts-numeric multiplication
#' follows the form of \code{Y*c(1,2)} if \code{Y} is a bivariate fts.
#' @return An object of class \code{\link{fts}}.
#'
#'
#' @param Y1 An object of class \code{\link{fts}} or numeric.
#' @param Y2 An object of class \code{\link{fts}} or numeric.
#' @examples
#' \dontrun{
#' require(Rfssa)
#' load_github_data("https://github.com/haghbinh/Rfssa/blob/main/data/Callcenter.RData")
#' D <- matrix(sqrt(Callcenter$calls), nrow = 240)
#' u <- seq(0, 1, length.out = 240) # Define domain of functional data
#' d <- 22 # number of basis elements
#' Y <- fts(list(D), list(list(d, "bspline")), list(u))
#' plot(Y)
#' Ytimes <- Y * Y # multiply the functional time series by itself
#' plot(Ytimes)
#' Ytimes2 <- Y * 2 # multiply every term in the fts by 2
#' plot(Ytimes2)
#' }
#'
#' @seealso \code{\link{fts}}
#' @useDynLib Rfssa
#' @export
"*.fts" <- function(Y1, Y2) {
  if (class(Y1) == "fts" && class(Y2) == "fts") {
    if (length(Y1@C) != length(Y2@C)) {
      stop("Functional time series are of different length.")
    }

    p <- length(Y1@C)
    grid <- Y1@grid
    basis <- Y1@B
    eval_fts <- list()
    new_grid <- list()
    for (j in 1:p) {
      if (max(max(Y1@B[[j]] - Y2@B[[j]])) >= 10^-14) {
        stop(paste("The basis matrix corresponding to variable ", as.character(j), " should be the same between fts objects.", sep = ""))
      }

      d <- ncol(basis[[j]])
      Y1_j_eval <- basis[[j]] %*% Y1@C[[j]]
      Y2_j_eval <- basis[[j]] %*% Y2@C[[j]]
      eval_fts[[j]] <- Y1_j_eval * Y2_j_eval
      N <- ncol(eval_fts[[j]])
      if (ncol(Y1@grid[[j]]) == 2) {
        x <- unique(Y1@grid[[j]][, 1])
        y <- unique(Y1@grid[[j]][, 2])
        new_fts_two_d <- array(data = NA, dim = c(length(x), length(y), N))
        for (n in 1:N) {
          count <- 1
          for (i_1 in 1:length(x)) {
            for (i_2 in 1:length(y)) {
              new_fts_two_d[i_1, i_2, n] <- eval_fts[[j]][count, n]
              count <- count + 1
            }
          }
        }

        eval_fts[[j]] <- new_fts_two_d
        new_grid[[j]] <- list(x, y)
      } else {
        new_grid[[j]] <- Y1@grid[[j]]
      }
    }
  } else if (class(Y1) == "fts" && class(Y2) == "numeric" || class(Y1) == "fts" && class(Y2) == "matrix") {
    Y2 <- as.numeric(Y2)
    p <- length(Y1@C)
    if (length(Y2) == 1) Y2 <- rep(Y2, p)
    grid <- Y1@grid
    basis <- Y1@B
    eval_fts <- list()
    new_grid <- list()
    for (j in 1:p) {
      d <- ncol(basis[[j]])
      Y1_j_eval <- basis[[j]] %*% Y1@C[[j]]
      eval_fts[[j]] <- Y1_j_eval * Y2[j]
      N <- ncol(eval_fts[[j]])
      if (ncol(Y1@grid[[j]]) == 2) {
        x <- unique(Y1@grid[[j]][, 1])
        y <- unique(Y1@grid[[j]][, 2])
        new_fts_two_d <- array(data = NA, dim = c(length(x), length(y), N))
        for (n in 1:N) {
          count <- 1
          for (i_1 in 1:length(x)) {
            for (i_2 in 1:length(y)) {
              new_fts_two_d[i_1, i_2, n] <- eval_fts[[j]][count, n]
              count <- count + 1
            }
          }
        }

        eval_fts[[j]] <- new_fts_two_d
        new_grid[[j]] <- list(x, y)
      } else {
        new_grid[[j]] <- Y1@grid[[j]]
      }
    }
  } else if (class(Y1) == "numeric" && class(Y2) == "fts" || class(Y1) == "matrix" && class(Y2) == "fts") {
    Y1 <- as.numeric(Y1)
    p <- length(Y2@C)
    if (length(Y1) == 1) Y1 <- rep(Y1, p)
    grid <- Y2@grid
    basis <- Y2@B
    eval_fts <- list()
    new_grid <- list()
    for (j in 1:p) {
      d <- ncol(basis[[j]])
      Y2_j_eval <- basis[[j]] %*% Y2@C[[j]]
      eval_fts[[j]] <- Y2_j_eval * Y1[j]
      N <- ncol(eval_fts[[j]])
      if (ncol(Y2@grid[[j]]) == 2) {
        x <- unique(Y2@grid[[j]][, 1])
        y <- unique(Y2@grid[[j]][, 2])
        new_fts_two_d <- array(data = NA, dim = c(length(x), length(y), N))
        for (n in 1:N) {
          count <- 1
          for (i_1 in 1:length(x)) {
            for (i_2 in 1:length(y)) {
              new_fts_two_d[i_1, i_2, n] <- eval_fts[[j]][count, n]
              count <- count + 1
            }
          }
        }

        eval_fts[[j]] <- new_fts_two_d
        new_grid[[j]] <- list(x, y)
      } else {
        new_grid[[j]] <- Y2@grid[[j]]
      }
    }
  }

  out <- Rfssa::fts(eval_fts, basis, new_grid)
  return(out)
}


#' Indexing into Functional Time Series
#'
#' An indexing method for functional time series (\code{\link{fts}}).
#' @return An object of class \code{\link{fts}}.
#'
#'
#' @param Y An object of class \code{\link{fts}}.
#' @param i The index value.
#' @examples
#' \dontrun{
#' require(Rfssa)
#' load_github_data("https://github.com/haghbinh/Rfssa/blob/main/data/Callcenter.RData")
#' D <- matrix(sqrt(Callcenter$calls), nrow = 240)
#' u <- seq(0, 1, length.out = 240) # Define domain of functional data
#' d <- 22 # number of basis elements
#' Y <- fts(list(D), list(list(d, "bspline")), list(u))
#' plot(Y)
#' plot(Y[10:15])
#' }
#'
#' @seealso \code{\link{fts}}
#' @useDynLib Rfssa
#' @export
"[.fts" <- function(Y, i = "index") {
  p <- length(Y@C)
  grid <- Y@grid
  basis <- Y@B
  eval_fts <- list()
  new_grid <- list()
  for (j in 1:p) {
    eval_fts[[j]] <- basis[[j]] %*% Y@C[[j]][, i]
    N <- ncol(eval_fts[[j]])

    if (ncol(Y@grid[[j]]) == 2) {
      x <- unique(Y@grid[[j]][, 1])
      y <- unique(Y@grid[[j]][, 2])
      new_fts_two_d <- array(data = NA, dim = c(length(x), length(y), N))
      for (n in 1:N) {
        count <- 1
        for (i_1 in 1:length(x)) {
          for (i_2 in 1:length(y)) {
            new_fts_two_d[i_1, i_2, n] <- eval_fts[[j]][count, n]
            count <- count + 1
          }
        }
      }

      eval_fts[[j]] <- new_fts_two_d
      new_grid[[j]] <- list(x, y)
    } else {
      new_grid[[j]] <- Y@grid[[j]]
    }
  }

  out <- Rfssa::fts(eval_fts, basis, new_grid)
  return(out)
}
