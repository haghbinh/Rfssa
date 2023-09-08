#' Weighted Correlation Matrix
#'
#' This function returns the weighted correlation (w-correlation) matrix for functional time series (\code{\link{funts}}) objects
#' that were reconstructed from functional singular spectrum analysis (\code{\link{fssa}}) objects.
#' @return A square matrix of w-correlation values for the reconstructed \code{\link{funts}} objects that were built from.
#' \code{\link{fssa}} components
#' @param U An object of class \code{\link{fssa}}.
#' @param groups A list or vector of indices which determines the grouping used for the reconstruction
#' in pairwise w-correlations matrix.
#' @examples
#' \dontrun{
#' ## Univariate FSSA Example on Callcenter data
#' ## Univariate functional singular spectrum analysis
#' data("Callcenter")
#' L <- 28
#' U <- fssa(Callcenter, L)
#' ufwcor <- fwcor(U = U, groups = list(1, 2, 3))
#' wplot(W = ufwcor)
#' }
#'
#' @seealso \code{\link{fssa}}, \code{\link{freconstruct}}, \code{\link{funts}}, \code{\link{wplot}}
#' @export
fwcor <- function(U, groups) {
  if (class(U[[1]])[[1]] != "list") out <- ufwcor(U, groups) else out <- mfwcor(U, groups)
  return(out)
}


# Univariate weighted correlation used to find weighted correlation matrix for the grouping stage of ufssa.

ufwcor <- function(U, groups) {
  if (is.numeric(groups)) groups <- as.list(groups)
  d <- length(groups)
  Q <- freconstruct(U, groups = groups)
  N <- U$N
  L <- U$L
  K <- N - L + 1L
  w <- 1L:N
  basis <- U$Y$B_mat[[1]]
  grid <- as.matrix(U$Y$argval[[1]])
  if (U$Y$dimSupp[[1]] == 1) {
    G <- onedG(A = basis, B = basis, grid = grid)
  } else {
    G <- twodG(A = basis, B = basis, grid = grid)
  }
  L1 <- min(L, K)
  K1 <- max(K, L)
  w[L1:K1] <- L1
  w[(K1 + 1L):N] <- N + 1L - ((K1 + 1L):N)
  out <- matrix(1L, nrow = d, ncol = d)
  for (i in 1L:(d - 1)) {
    for (j in (i + 1L):d) {
      out[i, j] <- winprod(Q[[i]]$coefs[[1]], Q[[j]]$coefs[[1]], w, G) / sqrt(winprod(Q[[i]]$coefs[[1]], Q[[i]]$coefs[[1]], w, G) * winprod(Q[[j]]$coefs[[1]], Q[[j]]$coefs[[1]], w, G))
    }
  }
  for (i in 2:d) for (j in 1:(i - 1)) out[i, j] <- out[j, i]
  return(out)
}

# Multivariate weighted correlation used to find weighted correlation matrix for the grouping stage of mfssa.

mfwcor <- function(U, groups) {
  if (is.numeric(groups)) groups <- as.list(groups)
  d <- length(groups)
  Q <- mfreconstruct(U, groups = groups)
  N <- U$N
  L <- U$L
  K <- N - L + 1L
  w <- 1L:N
  p <- length(U$Y$coefs)
  Y <- U$Y
  G <- list()
  for (i in 1:p) {
    grid <- as.matrix(Y$argval[[i]])
    if (Y$dimSupp[[i]] == 1) {
      G[[i]] <- t(onedG(A = Y$B_mat[[i]], B = Y$B_mat[[i]], grid = grid))
    } else {
      G[[i]] <- t(twodG(A = Y$B_mat[[i]], B = Y$B_mat[[i]], grid = grid))
    }
  }
  L1 <- min(L, K)
  K1 <- max(K, L)
  w[L1:K1] <- L1
  w[(K1 + 1L):N] <- N + 1L - ((K1 + 1L):N)
  wcor <- matrix(1L, nrow = d, ncol = d)
  for (i in 1L:(d - 1)) {
    Q_i <- Q[[i]]
    Q_i_l <- list()
    for (k in 1:p) {
      Q_i_l[[k]] <- Q_i$coefs[[k]]
    }
    for (j in (i + 1L):d) {
      Q_j <- Q[[j]]
      Q_j_l <- list()
      for (k in 1:p) {
        Q_j_l[[k]] <- Q_j$coefs[[k]]
      }
      wcor[i, j] <- mwinprod(Q_i_l, Q_j_l, w, G, p) / sqrt(mwinprod(Q_i_l, Q_i_l, w, G, p) * mwinprod(Q_j_l, Q_j_l, w, G, p))
    }
  }
  for (i in 2:d) for (j in 1:(i - 1)) wcor[i, j] <- wcor[j, i]
  return(wcor)
}
