# Utility functions

# This function Hankelize a d*K*L array.
fH <- function(C, d) {
  for (j in 1:d) {
    C[j, , ] <- H(C[j, , ])
  }
  return(C)
}

# Create Vectorize Lag vetors
# Get an d*N matrix, return an d*L matrix
lagvec_new <- function(coefs, L, i) {
  coefs[, i:(i + L - 1)]
}

# Projection of all lag vectors onto i-th functional eigen vector.
ufproj <- function(U, i, d) {
  L <- U$L
  K <- U$N - L + 1L
  Y <- U$Y
  B <- U$Y@B[[1]]
  u <- solve(t(B) %*% B) %*% t(B) %*% U[[i]]
  if (ncol(Y@grid[[1]]) == 1) {
    G <- onedG(A = B, B = B, grid = Y@grid[[1]])
  } else {
    G <- twodG(A = B, B = B, grid = Y@grid[[1]])
  }
  CX <- array(NA, dim = c(d, K, L))
  for (k in 1L:K) {
    x <- lagvec_new(Y@C[[1]], L, k)
    CX[, k, ] <- HLinprod(x, u, G) * u
  }
  return(CX)
}

# Right eigenvectors (return a K*d matrix)
# d is number of eigenfunctions
uV <- function(U, d) {
  L <- U$L
  N <- U$N
  K <- N - L + 1
  CX <- matrix(NA, nrow = K, ncol = d)
  basis <- U$Y@B[[1]]
  if (ncol(U$Y@grid[[1]]) == 1) {
    G <- onedG(A = basis, B = basis, grid = U$Y@grid[[1]])
  } else {
    G <- twodG(A = basis, B = basis, grid = U$Y@grid[[1]])
  }
  for (i in 1L:d) {
    u <- solve(t(basis) %*% basis) %*% t(basis) %*% U[[i]]
    lambd <- sqrt(U$values[i])
    for (k in 1L:K) {
      x <- lagvec_new(U$Y@C[[1]], L, k)
      CX[k, i] <- HLinprod(x, u, G) / lambd
    }
  }
  return(CX)
}
# right singular vectors for the multivariate case
mV <- function(U, d) {
  Y <- U$Y
  N <- U$N
  p <- length(U$Y@C)
  L <- U$L
  K <- N - L + 1L
  V <- matrix(
    nrow = K,
    ncol = d,
    data = NA
  )
  G <- list()
  # get inner product matrices
  for (i in 1L:p) {
    if (ncol(Y@grid[[i]]) == 1) {
      G[[i]] <- onedG(A = Y@B[[i]], B = Y@B[[i]], grid = Y@grid[[i]])
    } else {
      G[[i]] <- twodG(A = Y@B[[i]], B = Y@B[[i]], grid = Y@grid[[i]])
    }
  }
  for (i in 1L:d) {
    element <- U[[i]]
    u <- list()
    for (k in 1L:K) {
      x <- list()
      for (j in 1L:p) {
        u[[j]] <- solve(t(U$Y@B[[j]]) %*% U$Y@B[[j]]) %*% t(U$Y@B[[j]]) %*% element[[j]]
        x[[j]] <- lagvec_new(U$Y@C[[j]], L, k)
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
  K <- U$N - L + 1L
  Y <- U$Y
  u <- list()
  B <- U$Y@B
  p <- length(U$Y@C)
  for (j in 1:p) {
    u[[j]] <- solve(t(B[[j]]) %*% B[[j]]) %*% t(B[[j]]) %*% U[[i]][[j]]
  }
  C <- list()
  G <- list()
  for (j in 1:p) {
    d <- ncol(Y@B[[j]])
    C[[j]] <- array(NA, dim = c(d, K, L))
    if (ncol(Y@grid[[j]]) == 1) {
      G[[j]] <- t(onedG(A = Y@B[[j]], B = Y@B[[j]], grid = Y@grid[[j]]))
    } else {
      G[[j]] <- t(twodG(A = Y@B[[j]], B = Y@B[[j]], grid = Y@grid[[j]]))
    }
  }
  # define HpL lag vector
  for (k in 1:K) {
    x <- list()
    for (j in 1:p) {
      x[[j]] <- lagvec_new(Y@C[[j]], L, k)
    }
    # build operator
    for (j in 1:p) {
      C_jx <- C[[j]]
      C_jx[, k, ] <- HpLinprod(x, u, G, p) * u[[j]]
      C[[j]] <- C_jx
    }
  }
  return(C)
}

onedG <- function(A, B, grid, method = "trapezoidal") {
  M <- nrow(grid)
  dx <- grid[2] - grid[1]
  if (method == "rectangular") {
    G <- dx * t(B) %*% A
  } else {
    D_diag <- c(1, rep(x = 2, (M - 2)), 1)
    D <- matrix(data = 0, nrow = M, ncol = M)
    diag(D) <- D_diag
    G <- (dx / 2) * t(B) %*% D %*% A
  }
  return(G)
}

twodG <- function(A, B, grid, method = "trapezoidal") {
  if (grid[2, 1] == grid[1, 1]) { # case for by columns
    dy <- grid[2, 2] - grid[1, 2]
    M_1 <- length(which(grid[, 1] == grid[1, 1]))
    dx <- grid[(M_1 + 1), 1] - grid[1, 1]
    M_2 <- length(which(grid[, 2] == grid[1, 2]))
  } else { # case for by rows
    dx <- grid[2, 1] - grid[1, 1]
    M_2 <- length(which(grid[, 2] == grid[1, 2]))
    dy <- grid[(M_2 + 1), 2] - grid[1, 2]
    M_1 <- length(which(grid[, 1] == grid[1, 1]))
  }
  M <- M_1 * M_2
  if (method == "rectangular") {
    G <- (dx * dy) * t(B) %*% A
  } else {
    a_1 <- c(1, rep(x = 2, (M_1 - 2)), 1)
    a_2 <- c(2, rep(x = 4, (M_1 - 2)), 2)
    D_diag <- c(a_1, rep(a_2, (M_2 - 2)), a_1)
    D <- matrix(data = 0, nrow = M, ncol = M)
    diag(D) <- D_diag
    G <- (dx * dy / 4) * t(B) %*% D %*% A
  }

  return(G)
}



init_basis_check <- function(B) {
  if (!is.basis(B) & is.list(B)) {
    if (!all(sapply(B, is.basis)) | !all(sapply(B, is.basis))) {
      stop("All elements of the basis list must be `basisfd` objects.")
    }
  }
  if(!is.basis(basis) & !is.list(basis)){
    stop("The basis must be a `basisfd` object or list of `basisfd` objects.")
  }
}
