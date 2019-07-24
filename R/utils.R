# Code built by Jordan Trinka and Mehdi Maadooliat of Marquette University and Hossein Haghbin of Persian Gulf University



# This funcrtion Henkelize an d*K*L array.
fH <- function(C, d) {
  library("Rcpp")
  sourceCpp("H.cpp")
  for (j in 1:d) C[j, , ] <- H(C[j, , ])
  return(C)
}

# Create Vectorize Lag vetors
# Get an d*N matrix, return an d*L matrix
lagvec_new <- function(coefs, L, i)  coefs[, i:(i + L - 1)]

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
# right singular vecjtors for the multivariate case
mV<-function(U,values,K,Y,A,p,L){
  sourceCpp("HpLinprod.cpp")
  V <- matrix(nrow=K,ncol=length(values),data=0)
  for(i in 1:length(values)){
    element <- U[[i]]
    u <- list()
    for(k in 1:K){
      x <- list()
      for(j in 1:p){
        u[[j]] <- element[[j]]$coefs
        x[[j]] <- lagvec_new(Y[[j]]$coefs,L,k)
      }
      V[k,i] <- HpLinprod(u,x,A,p)/sqrt(values[i])
    }
  }
  return(V)
}

# Rank one approximajtion
mfproj <- function(out, i, K, L, Y){
  # defining pieces of multivariate rank one approximation
  sourceCpp("HpLinprod.cpp")
  u <- list()
  pc_1 <- out$U[[1]]
  p <- length(pc_1)
  for(j in 1 : p){
    element<-out$U[[i]]
    u[[j]]=element[[j]]$coefs
  }
  C = list()
  A = list()
  for(j in 1:p){

    d=nrow(pc_1[[j]]$coefs)
    C[[j]]=array(NA, dim = c(d, K, L))
    A[[j]] = inprod(Y[[j]]$basis,Y[[j]]$basis)

  }
  # define HpL lag vector
  for(k in 1:K){
    x <- list()
    for(j in 1:p){
      x[[j]]=lagvec_new(Y[[j]]$coefs,L,k)

    }
    #build operajtor
    for(j in 1:p){
      C_jx <- C[[j]]
      C_jx[, k, ] <- HpLinprod(x,u,A,p)*u[[j]]
      C[[j]] <- C_jx

    }

  }
  return(C)
}
