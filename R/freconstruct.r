#' Reconstruction Stage of Functional Singular Spectrum Analysis
#'
#' Reconstruct univariate or multivariate functional time series
#' (\code{\link{funts}}) objects from functional singular spectrum analysis
#' (\code{\link{fssa}}) objects, including Grouping and Hankelization steps.
#' This function performs the reconstruction step for either univariate
#' functional singular spectrum analysis (ufssa) or multivariate
#' functional singular spectrum analysis (mfssa), depending on the input.
#'
#' @param U an object of class \code{\link{fssa}}.
#' @param groups a list of numeric vectors, each vector includes indices
#' of elementary components of a group used for reconstruction.
#'
#' @return A named list of objects of class \code{\link{funts}} that are
#' reconstructed according to the specified groups and a numeric vector
#' of eigenvalues.
#'
#' @note Refer to \code{\link{fssa}} for an example on how to run this
#' function starting from \code{\link{fssa}} objects.
#'
#' @seealso \code{\link{fssa}}, \code{\link{funts}}
#'
#' @export
freconstruct <- function(U, groups = as.list(1L:10L)) {
  if (class(U[[1]])[[1]] != "list") out <- ufreconstruct(U, groups) else out <- mfreconstruct(U, groups)
  return(out)
}




#--------------------------------------ufreconstruct----------------------------------------------------


# Reconstruction stage (including Hankelization) of univariate functional singular spectrum analysis

ufreconstruct <- function(U, groups = as.list(1L:10L)) {
  Y <- U$Y
  N <- Y$N
  basisobj <- Y$basis[[1]]
  time_st <- Y$time[1]
  time_en <- Y$time[N]
  basis <- Y$B_mat[[1]]
  d <- ncol(Y$B_mat[[1]])
  L <- U$L
  K <- N - L + 1L
  m <- length(groups)
  out <- list()
  for (i in 1L:m) {
    Cx <- matrix(NA, nrow = d, ncol = N)
    g <- groups[[i]]
    S <- 0L
    for (j in 1L:length(g)) S <- S + ufproj(U, g[j], d)
    S <- fH(S, d)
    Cx[, 1L:L] <- S[, 1L, ]
    Cx[, L:N] <- S[, , L]
    recon_out <- basis %*% Cx
    if (Y$dimSupp[[1]] == 2) {
      x <- unique(Y$argval[[1]][, 1])
      y <- unique(Y$argval[[1]][, 2])
      recon_two_d <- array(data = NA, dim = c(length(x), length(y), N))
      for (n in 1:N) {
        count <- 1
        for (i_1 in 1:length(x)) {
          for (i_2 in 1:length(y)) {
            recon_two_d[i_1, i_2, n] <- recon_out[count, n]
            count <- count + 1L
          }
        }
      }
      recon_out <- recon_two_d
      new_grid <- list(x, y)
    } else {
      new_grid <- Y$argval[[1]]
    }
    funts_out <- funts(X = recon_out, basisobj = basisobj, argval = new_grid, start = time_st, end = time_en)
    out[[i]] <- funts_out
  }
  out$values <- sqrt(U$values)
  return(out)
}


#--------------------------------------mfreconstruct----------------------------------------------------

# Reconstruction stage (including Hankelization) of multivariate functional singular spectrum analysis.
mfreconstruct <- function(U, groups = as.list(1L:10L)) {
  Y <- U$Y
  N <- Y$N
  basisobj <- Y$basis
  time_st <- Y$time[1]
  time_en <- Y$time[N]
  L <- U$L
  K <- N - L + 1
  basis <- Y$B_mat
  m <- length(groups)
  recon_out <- list()
  new_grid <- list()
  p <- length(Y$dimSupp)
  # Loop over groups
  for (i in 1L:m) {
    recon_group <- list()
    C <- list()
    S <- list()
    g <- groups[[i]]
    # build reconstructions
    for (j in 1:p) {
      d <- ncol(basis[[j]])
      C[[j]] <- matrix(NA, nrow = d, ncol = N)
      S[[j]] <- 0L
      for(k in 1L:length(g)){
        S[[j]] <- S[[j]] + mfproj(U, g[k])[[j]]
      }
      S[[j]] <- fH(S[[j]], d)
      C_jx <- C[[j]]
      S_jx <- S[[j]]
      C_jx[, 1L:L] <- S_jx[, 1L, ]
      C_jx[, L:N] <- S_jx[, , L]
      recon_group[[j]] <- basis[[j]] %*% C_jx
      if (Y$dimSupp[[j]] == 2) {
        x <- unique(Y$argval[[j]][, 1])
        y <- unique(Y$argval[[j]][, 2])
        recon_two_d <- array(data = NA, dim = c(length(x), length(y), N))
        for (n in 1:N) {
          count <- 1
          for (i_1 in 1:length(x)) {
            for (i_2 in 1:length(y)) {
              recon_two_d[i_1, i_2, n] <- recon_group[[j]][count, n]
              count <- count + 1
            }
          }
        }

        recon_group[[j]] <- recon_two_d
        new_grid[[j]] <- list(x, y)
      } else {
        new_grid[[j]] <- Y$argval[[j]]
      }
    }

    # output the reconstructions
    funts_out <- funts(X = recon_group, basisobj = basisobj, argval = new_grid, start = time_st, end = time_en)
    recon_out[[i]] <- funts_out
  }
  recon_out$values <- U$values
  return(recon_out)
}

