# Utility functions

# This function Hankelize a d*K*L array.
fH <- function(C, d) {
  for (j in 1:d)
    C[j, ,] <- H(C[j, ,])
  return(C)
}

# Create Vectorize Lag vetors
# Get an d*N matrix, return an d*L matrix
lagvec_new <- function(coefs, L, i)
  coefs[, i:(i + L - 1)]

# Projection of all lag vectors onto i-th functional eigen vector.
ufproj <- function(U, i, d) {
  L <- U$L
  K <- U$N - L +1L
  Y <- U$Y
  u <- U[[i]]$coefs
  basis <- Y[[1]]$basis
  G <- inprod(basis, basis)
  CX <- array(NA, dim = c(d, K, L))
  for (k in 1L:K) {
    x <- lagvec_new(Y[[1]]$coefs, L, k)
    CX[, k,] <- HLinprod(x, u, G) * u
  }
  return(CX)
}

# Right eigenvectors (return a K*d matrix)
# d is number of eigenfunctions
uV <- function(U,d) {
  L <- U$L
  N <- U$N
  K <- N-L+1
  CX <- matrix(NA, nrow = K, ncol = d)
  basis <- U$Y[[1]]$basis
  G <- inprod(basis, basis)
  for (i in 1L:d) {
    u <- U[[i]]$coefs
    lambd <- sqrt(U$values[i])
    for (k in 1L:K) {
      x <- lagvec_new(U$Y[[1]]$coefs, L, k)
      CX[k, i] <- HLinprod(x, u, G) / lambd
    }
  }
  return(CX)
}
# right singular vectors for the multivariate case
mV <- function(U, d) {
  Y <- U$Y
  N <- U$N
  p <- Y$p
  L <- U$L
  K <- N - L + 1L
  V <- matrix(nrow = K,
              ncol = d,
              data = NA)
  G <- list()
  # get inner product matrices
  for (i in 1L:p) {
    basis <- U$Y[[i]]$basis
    G[[i]] <- inprod(basis,basis)
  }
  for (i in 1L:d) {
    element <- U[[i]]
    u <- list()
    for (k in 1L:K) {
      x <- list()
      for (j in 1L:p) {
        u[[j]]  <- element[[j]]$coefs
        x[[j]] <- lagvec_new(Y[[j]]$coefs, L, k)
      }
      V[k, i] <- HpLinprod(u, x, G, p) / sqrt(U$values[i])
    }
  }
  return(V)
}

# Rank one approximation
mfproj <- function(U, i) {
  # defining pieces of multivariate rank one approximation
  L <- U$L
  K <- U$N-L+1L
  Y <- U$Y
  u <- list()
  pc_1 <- U[[1]]
  p <- length(pc_1)
  for (j in 1:p) {
    element <- U[[i]]
    u[[j]] <- element[[j]]$coefs
  }
  C <- list()
  G <- list()
  for (j in 1:p) {
    d <- nrow(pc_1[[j]]$coefs)
    C[[j]] <- array(NA, dim = c(d, K, L))
    G[[j]] <- inprod(Y[[j]]$basis, Y[[j]]$basis)
  }
  # define HpL lag vector
  for (k in 1:K) {
    x <- list()
    for (j in 1:p) {
      x[[j]] <- lagvec_new(Y[[j]]$coefs, L, k)
    }
    #build operator
    for (j in 1:p) {
      C_jx <- C[[j]]
      C_jx[, k,] <- HpLinprod(x, u, G, p) * u[[j]]
      C[[j]] <- C_jx
    }
  }
  return(C)
}
