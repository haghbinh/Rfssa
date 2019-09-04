# Univariate and multivariate weighted correlation used to find weighted correlation matrix for the grouping stage
# of ufssa and mfssa.

ufwcor <- function(U, d) {
  Q <- freconstruct(U, group = as.list(1:d))
  N <- U$N
  L <- U$L
  K <- N - L + 1L
  w <- 1L:N
  basis <- U$Y[[1]]$basis
  G <- inprod(basis, basis)
  L1 <- min(L, K)
  K1 <- max(K, L)
  w[L1:K1] <- L1
  w[(K1 + 1L):N] <- N + 1L - ( (K1 +1L):N)
  out <- matrix(1L, nrow = d, ncol = d)
  for (i in 1L:(d - 1)) {
    for (j in (i + 1L):d) {
      out[i, j] <- winprod(Q[[i]][[1]]$coefs, Q[[j]][[1]]$coefs, w, G)/sqrt(winprod(Q[[i]][[1]]$coefs,Q[[i]][[1]]$coefs, w, G) * winprod(Q[[j]][[1]]$coefs,Q[[j]][[1]]$coefs, w, G))
    }
  }
  for (i in 2:d) for (j in 1:(i - 1)) out[i, j] <- out[j, i]
  return(out)
}

mfwcor <- function(U, d) {
  Q <- mfreconstruct(U, group = as.list(1:d))
  N <- U$N
  L <- U$L
  K <- N - L + 1L
  w <- 1L:N
  p = U$Y$p
  Y = U$Y
  G = list()
  for(j in 1:p){
    G[[j]] = inprod(Y[[j]]$basis,Y[[j]]$basis)
  }
  L1 <- min(L, K)
  K1 <- max(K, L)
  w[L1:K1] <- L1
  w[(K1 + 1L):N] <- N + 1L - ( (K1 +1L):N)
  wcor <- matrix(1L, nrow = d, ncol = d)
  for (i in 1L:(d - 1)){
    Q_i=Q[[i]]
    Q_i_l <- list()
    for(k in 1:p){
      Q_i_l[[k]]=Q_i[[k]]$coefs
    }
    for (j in (i + 1L):d){
      Q_j=Q[[j]]
      Q_j_l <- list()
      for(k in 1:p){
        Q_j_l[[k]]=Q_j[[k]]$coefs
      }
      wcor[i, j] <- mwinprod(Q_i_l, Q_j_l, w, G, p)/sqrt(mwinprod(Q_i_l,Q_i_l, w, G, p) * mwinprod(Q_j_l,Q_j_l, w, G, p))
    }
  }
  for (i in 2:d) for (j in 1:(i - 1)) wcor[i, j] <- wcor[j, i]
  return(wcor)
}

#' W-correlation matrix
#'
#' This function evaluate the Wcorrelation plot for fssa
#' @return A squared matrix of W-correlation values.
#' @param U in the input is an object of class \code{fssa}.
#' @param d is the number of elementary components
#' in pairwise W-correlations matrix.
#' @seealso \code{\link{fssa}}
#' @export
fwcor <- function(U, d) {
  if(is.fd(U[[1]])) out <- ufwcor(U, d) else out <- mfwcor(U, d)
  return(out)
}
