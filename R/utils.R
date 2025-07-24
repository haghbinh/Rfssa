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

#Check the provided basis when creating funts object
check_basis <- function(basis){
  G = t(basis)%*%basis
  lambda_min = min(eigen(x = G,symmetric = TRUE)$values)
  return(lambda_min)
}

# Projection of all lag vectors onto i-th functional eigen vector.
ufproj <- function(U, i, d) {
  L <- U$L
  Y <- U$Y
  N <- Y$N
  K <- N - L + 1L
  B <- Y$B_mat[[1]]
  u <- solve(t(B) %*% B) %*% t(B) %*% U[[i]]
  grid <- as.matrix(Y$argval[[1]])
  basis <- Y$B_mat[[1]]
  if (Y$dimSupp[[1]] == 1) {
    G <- onedG(A = basis, B = basis, grid = grid)
  } else {
    G <- twodG(A = basis, B = basis, grid = grid)
  }
  CX <- array(NA, dim = c(d, K, L))
  for (k in 1L:K) {
    x <- lagvec_new(Y$coefs[[1]], L, k)
    CX[, k, ] <- HLinprod(x, u, G) * u
  }
  return(CX)
}

# Right eigenvectors (return a K*d matrix)
# d is number of eigenfunctions
uV <- function(U, d) {
  L <- U$L
  Y <- U$Y
  N <- Y$N
  K <- N - L + 1L
  CX <- matrix(NA, nrow = K, ncol = d)
  basis <- Y$B_mat[[1]]
  grid <- as.matrix(Y$argval[[1]])
  if (Y$dimSupp[[1]] == 1) {
    G <- onedG(A = basis, B = basis, grid = grid)
  } else {
    G <- twodG(A = basis, B = basis, grid = grid)
  }
  for (i in 1L:d) {
    u <- solve(t(basis) %*% basis) %*% t(basis) %*% U[[i]]
    lambd <- sqrt(U$values[i])
    for (k in 1L:K) {
      x <- lagvec_new(U$Y$coefs[[1]], L, k)
      CX[k, i] <- HLinprod(x, u, G) / lambd
    }
  }
  return(CX)
}
# right singular vectors for the multivariate case
mV <- function(U, d) {
  Y <- U$Y
  N <- Y$N
  p <- length(Y$dimSupp)
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
    grid <- as.matrix(Y$argval[[i]])
    basis <- Y$B_mat[[i]]
    if (Y$dimSupp[[i]] == 1) {
      G[[i]] <- onedG(A = basis, B = basis, grid = grid)
    } else {
      G[[i]] <- twodG(A = basis, B = basis, grid = grid)
    }
  }
  for (i in 1L:d) {
    element <- U[[i]]
    u <- list()
    for (k in 1L:K) {
      x <- list()
      for (j in 1L:p) {
        u[[j]] <- solve(t(U$Y$B_mat[[j]]) %*% U$Y$B_mat[[j]]) %*% t(U$Y$B_mat[[j]]) %*% element[[j]]
        x[[j]] <- lagvec_new(U$Y$coefs[[j]], L, k)
      }
      V[k, i] <- HpLinprod(u, x, G, p) / sqrt(U$values[i])
    }
  }
  return(V)
}

# Rank one approximation
mfproj <- function(U, i) {
  # defining pieces of multivariate rank one approximation
  Y <- U$Y
  N <- Y$N
  p <- length(Y$dimSupp)
  L <- U$L
  K <- N - L + 1L
  B <- Y$B_mat
  u <- list()
  for (j in 1:p) {
    u[[j]] <- solve(t(B[[j]]) %*% B[[j]]) %*% t(B[[j]]) %*% U[[i]][[j]]
  }
  C <- list()
  G <- list()
  for (j in 1:p) {
    d <- ncol(B[[j]])
    C[[j]] <- array(NA, dim = c(d, K, L))
    grid <- as.matrix(Y$argval[[j]])
    basis <- B[[j]]
    if (Y$dimSupp[[j]] == 1) {
      G[[j]] <- onedG(A = basis, B = basis, grid = grid)
    } else {
      G[[j]] <- twodG(A = basis, B = basis, grid = grid)
    }
  }
  # define HpL lag vector
  for (k in 1:K) {
    x <- list()
    for (j in 1:p) {
      x[[j]] <- lagvec_new(Y$coefs[[j]], L, k)
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


# A very Simple Function that evaluates Empirical Basis in any grid points by linear approximation.
eval.empb <- function(evalarg, basisobj) {
  if (!is.null(attr(basisobj, "grids"))) {
    grids <- attr(basisobj, "grids")
  } else if (!is.null(attr(basisobj, "rangeval"))) {
    r <- attr(basisobj, "rangeval")
    grids <- seq(r[1], r[2], len = nrow(basisobj))
  } else {
    grids <- seq(0, 1, len = nrow(basisobj))
  }
  fun.b <- apply(basisobj, 2, approxfun, x = grids)
  b <- NULL
  for (i in 1:length(fun.b)) b <- cbind(b, fun.b[[i]](evalarg))
  attr(b, "rangeval") <- range(evalarg)
  attr(Bs1, "grids") <- evalarg
  return(b)
}

# init_basis_check <- function(B) {
#   if (!is.basis(B) & is.list(B)) {
#     if (!all(sapply(B, is.basis)) | !all(sapply(B, is.basis))) {
#       stop("All elements of the basis list must be `basisfd` objects.")
#     }
#   }
#   if (!is.basis(B) & !is.list(B)) {
#     stop("The basis must be a `basisfd` object or list of `basisfd` objects.")
#   }
# }
