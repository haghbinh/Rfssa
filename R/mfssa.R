# Code built by Jordan Trinka and Mehdi Maadooliat of Marquette University and Hossein Haghbin of Persian Gulf University
# MFSSA Decomposition
mfssa <- function(Y, L = floor(dim(Y$coefs)[2L]/2L)){
  # get c plus plus code
  p <- Y$p
  d <- matrix(c(0,Y$d),nrow = 1L, ncol = (p+1L))
  N <- Y$N
  B <- list()
  A <- list()
  # get inner producjt matrices
  for(i in 1:p){
    B[[i]] <- inprod(Y[[i]],Y[[i]]$basis)
    A[[i]] <- inprod(Y[[i]]$basis,Y[[i]]$basis)
  }
  # Find the proper inner product matrices for j_k variables
  d_tilde <- sum(d)/L
  K <- N - L + 1L
  # change from ranges to shifter
  shifter <- matrix(nrow = 2, ncol = (p+1L), data=0L)
  shifter[,2L] <- c(1L,d[2L])
  if(p > 1L){
    for(i in 2L:p){
      shifter[1L,i+1L]=shifter[2L,i]+1L
      shifter[2L,i+1L]=shifter[2L,i]+d[i+1L]
    }
  }
  # find the desired majtrices
  S_0 <- SSM(K, L, d_tilde, p, B, shifter)
  G <- Gramm(K,L,p,d_tilde,A,shifter,d)
  S <- solve(G)%*%S_0 # S matrix which parameterizes var/cov op.
  Q <- eigen(S)
  coefs <- Re(Q$vectors)
  p_c <- list()
  r <- sum(Re(Q$values) > 0.001)
  values <- Re(Q$values[1L:r])
  out <- list()
  for(i in 1L:(r)){
    my_pcs <- list(NA)
    for(j in 1L: p){
      my_pcs[[j]] <- fd(Cofmat((d[j+1L]/L), L, coefs[(shifter[1L,(j+1L)]:shifter[2L,(j+1L)]),i]),Y[[j]]$basis)
    }
    out[[i]] <- my_pcs
  }

  out$values <- values
  out$L <- L
  out$N <- N
  out$Y <- Y
  class(out)="mfssa"
  return(out)
}
