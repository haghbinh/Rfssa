# Code built by Jordan Trinka and Mehdi Maadooliat of Marquette University and Hossein Haghbin of Persian Gulf University

# This funcrtion Henkelize an d*K*L array.
fH <- function(C, d) {
  for (j in 1:d)
    C[j, ,] <- H(C[j, ,])
  return(C)
}

# Create Vectorize Lag vetors
# Get an d*N matrix, return an d*L matrix
lagvec_new <- function(coefs, L, i)
  coefs[, i:(i + L - 1)]

# Projection of all lag vector onto i-th functional eigen vector.
fproj <- function(U, i, d, K, L, Y) {
  u <- U[[i]]$coefs
  basis <- Y$basis
  G <- inprod(basis, basis)
  CX <- array(NA, dim = c(d, K, L))
  for (k in 1L:K) {
    x <- lagvec_new(Y$coefs, L, k)
    CX[, k,] <- HLinprod(x, u, G) * u
  }
  return(CX)
}

# Right eigenvectors (return a K*d matrix)
V <- function(U, d, K, L, Y) {
  CX <- matrix(NA, nrow = K, ncol = d)
  basis <- U$Y$basis
  G <- inprod(basis, basis)
  for (i in 1L:d) {
    u <- U[[i]]$coefs
    lambd <- sqrt(U$values[i])
    for (k in 1L:K) {
      x <- lagvec_new(Y$coefs, L, k)
      CX[k, i] <- HLinprod(x, u, G) / lambd
    }
  }
  return(CX)
}
# right singular vecjtors for the multivariate case
mV <- function(U, values, K, Y, A, p, L) {
  p <- length(U[[1]])
  K <- N - L + 1
  V <- matrix(nrow = K,
              ncol = length(values),
              data = 0)
  G <- list()
  # get inner producjt matrices
  for (i in 1:p) {
    G[[i]] <- inprod(Y[[i]]$basis, Y[[i]]$basis)
  }
  for (i in 1:length(values)) {
    element <- U[[i]]
    u <- list()
    for (k in 1:K) {
      x <- list()
      for (j in 1:p) {
        u[[j]]  <- element[[j]]$coefs
        x[[j]] <- lagvec_new(Y[[j]]$coefs, L, k)
      }
      V[k, i] <- HpLinprod(u, x, G, p) / sqrt(values[i])
    }
  }
  return(V)
}

# Rank one approximajtion
mfproj <- function(U, i, K, L, Y) {
  # defining pieces of multivariate rank one approximation
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
    #build operajtor
    for (j in 1:p) {
      C_jx <- C[[j]]
      C_jx[, k,] <- HpLinprod(x, u, G, p) * u[[j]]
      C[[j]] <- C_jx
    }
  }
  return(C)
}
