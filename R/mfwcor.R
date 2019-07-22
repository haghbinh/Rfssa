mfwcor <- function(out, d) {
  source("mfreconstruct.R")
  library("Rcpp")
  sourceCpp("mwinprod.cpp")
  Q <- mfreconstruct(out, group = as.list(1:d))
  N <- out$N
  L <- out$L
  K <- N - L + 1L
  w <- 1L:N
  p = length(out$U[[1]])
  Y = out$Y
  A = list()
  for(j in 1:p){
    A[[j]] = inprod(Y[[j]]$basis,Y[[j]]$basis)
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
      wcor[i, j] <- mwinprod(Q_i_l, Q_j_l, w, A, p)/sqrt(mwinprod(Q_i_l,Q_i_l, w, A, p) * mwinprod(Q_j_l,Q_j_l, w, A, p))
    }
  }
  for (i in 2:d) for (j in 1:(i - 1)) wcor[i, j] <- wcor[j, i]
  return(wcor)
}