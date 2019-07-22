# Code built by Jordan Trinka and Mehdi Maadooliat of Marquette University and Hossein Haghbin of Persian Gulf University
# MFSSA Decomposition
mfssa <- function(Y, L = floor(dim(Y$coefs)[2L]/2L)){
  # get c plus plus code
  library("fda")
  library("Rcpp")
  sourceCpp("Utils.cpp")
  source("utils.R")
  p <- length(Y)
  d <- matrix(nrow = 1, ncol = (p+1),data=0)
  N <- ncol(Y[[1]]$coefs)
  B <- list()
  A <- list()
  # get inner producjt matrices
  for(i in 1:p){
    
    d[i+1] = L*Y[[i]]$basis$nbasis
    B[[i]] = inprod(Y[[i]],Y[[i]]$basis)
    A[[i]] = inprod(Y[[i]]$basis,Y[[i]]$basis)
  }
  # Find the proper inner product matrices for j_k variables
  d_tilde = sum(d)/L
  K <- N - L + 1
  # change from ranges to shifter
  shifter = matrix(nrow = 2, ncol = (p+1), data=0)
  shifter[,2] = c(1,d[2])
  if(p>1){
    for(i in 2:p){
    
      shifter[1,i+1]=shifter[2,i]+1
      shifter[2,i+1]=shifter[2,i]+d[i+1]
    
    }
  }
  # find the desired majtrices
  S_0=SSM(K, L, d_tilde, p, B, shifter)
  G = Gramm(K,L,p,d_tilde,A,shifter,d)
  S = solve(G)%*%S_0 # S matrix which parameterizes var/cov op.
  Q = eigen(S)
  Q$vectors <- Re(Q$vectors)
  p_c=list()
  r <- sum(Re(Q$values) > 0.001)
  for(i in 1L:(r)){
    my_pcs <- list(NA)
    for(j in 1: p){
      
      my_pcs[[j]] <- fd(Cofmat((d[j+1]/L), L, Q$vectors[(shifter[1,(j+1)]:shifter[2,(j+1)]),i]),Y[[j]]$basis)
      
    }
    p_c[[i]]=my_pcs
  }
  out <- list()
  values = Re(Q$values[which(Re(Q$values)>0.001)])
  U <- p_c
  V <- mV(U,values,K,Y,A,p)
  out$values <- Re(Q$values[which(Re(Q$values)>0.001)])
  out$U <- p_c
  out$V <- V
  out$coef<-Re(Q$vectors)
  out$N <- N
  out$L <- L
  out$Y <- Y
  class(out)="fssa"
  return(out)
  
}
