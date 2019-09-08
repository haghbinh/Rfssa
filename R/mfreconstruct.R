# Reconstruction stage (including Hankelization) of multivariate functional singular spectrum analysis.
mfreconstruct <-function(U, group = as.list(1L:10L)){
  N <- U$N
  Y <- U$Y
  p <- U$Y$p
  L <- U$L
  K <- N-L+1
  basis <- list()
  m <- length(group)
  for(j in 1: p){
    basis[[j]] <- U$Y[[j]]$basis
  }
  recon_out <- list()
  # Loop over groups
  for(i in 1L:m){
    recon_out_j <- list()
    C <- list()
    S <- list()
    element <- U[[1]]
    for(j in 1:p){
      d <- nrow(element[[j]]$coefs)
      C[[j]] <- matrix(NA, nrow = d, ncol = N)
      S[[j]] <- 0L
    }
    g <- group[[i]]
    for(k in 1L:length(g)){
      # projtect onto the space spanned by each p.c.
      projection <- mfproj(U,g[k])
      for(j in 1:p){
        S[[j]] <- S[[j]]+projection[[j]]
      }
    }
    # build reconstructions
    for(j in 1:p){
      d <- nrow(element[[j]]$coefs)
      S[[j]] <- fH(S[[j]],d)
      C_jx <- C[[j]]
      S_jx <- S[[j]]
      C_jx[,1L:L] <- S_jx[, 1L,]
      C_jx[,L:N] <- S_jx[,,L]
      recon_out_j[[j]] <- fd(C_jx, basis[[j]])
    }
    # output the reconstructions
    recon_out[[i]] <- fts(recon_out_j,time=Y$time)
  }
  recon_out$values <- U$values
  return(recon_out)
}

