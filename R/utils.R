# Projection of all lag vector onto i-th functional eigen vector.
fproj <- function(U, i, d, K, L, Y) {
    u <- U[[i]]$coefs
    basis <- U$Y$basis
    G <- inprod(basis,basis)
    CX <- array(NA, dim = c(d, K, L))
    for (k in 1L:K) {
        x <- lagvec_new(Y$coefs, L, k)
        CX[, k, ] <- HLinprod(x, u,G) * u
    }
    return(CX)
    }

# Right eigenvectors (return a K*d matrix)
V <- function(U, d, K, L, Y) {
  CX <- matrix(NA, nrow = K, ncol = d)
  basis <- U$Y$basis
  G <- inprod(basis,basis)
  for(i in 1L:d){
    u <- U[[i]]$coefs
    lambd <- sqrt(U$values[i])
    for (k in 1L:K) {
      x <- lagvec_new(Y$coefs, L, k)
      CX[k,i] <- HLinprod(x, u,G)/lambd
    }
  }
  return(CX)
}
# This funcrtion Henkelize an d*K*L array.
fH <- function(C, d) {
    for (j in 1L:d) C[j, , ] <- H(C[j, , ])
    return(C)
}

# Create Vectorize Lag vetors
# Get an d*N matrix, return an d*L matrix
lagvec_new <- function(coefs, L, i)  coefs[, i:(i + L - 1L)]


