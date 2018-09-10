winprod <- function(x, y, w, G) {
    Xcoef <- x$coefs
    Ycoef <- y$coefs
    N <- ncol(x$coefs)
    s <- NA
    for (i in 1:N) s[i] <- t(Xcoef[, i]) %*% G %*% Ycoef[, i]
    return(sum(s))
}

# Projection of all lag vector onto i-th functional eigen vector.
fproj <- function(U, i, d, K, L,
    Y) {
    u <- U[[i]]$coefs
    CX <- array(NA, dim = c(d, K, L))
    for (k in 1:K) {
        x <- lagvec_new(Y$coefs, L, k)
        CX[, k, ] <- HLinprod(x, u) * u
    }
    return(CX)
    }

# This funcrtion Henkelize an d*K*L array.
fH <- function(C, d) {
    for (j in 1:d) C[j, , ] <- H(C[j, , ])
    return(C)
}

# Create Vectorize Lag vetors
# Get an d*N matrix, return an d*L matrix
lagvec_new <- function(coefs, L, i)  coefs[, i:(i + L - 1)]

# Inner product between two elements of H^L.
HLinprod <- function(x,y){
  L <- ncol(x)
  s <- 0
  for(j in 1:L) s <- s+as.numeric(t(x[,j])%*%y[,j])
  return(s)
}

