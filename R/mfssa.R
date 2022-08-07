# Embedding and decomposition stages of multivariate functional singular spectrum analysis.
mfssa <- function(Y, L, ntriples) {
  # get c plus plus code
  p <- length(Y@C)
  N <- ncol(Y@C[[1]])
  C_tilde <- list()
  G_1 <- list()
  shifter <- matrix(nrow = 2, ncol = (p + 1L), data = 0L)
  Y_d <- matrix(data = 0, nrow = 1, ncol = p)
  # get inner product matrices
  for (j in 1:p) {

    if (ncol(Y@grid[[j]]) == 1) {
      C_tilde[[j]] <- t(onedG(A = Y@B[[j]] %*% Y@C[[j]], B = Y@B[[j]], grid = Y@grid[[j]])) # old B
      G_1[[j]] <- t(onedG(A = Y@B[[j]], B = Y@B[[j]], grid = Y@grid[[j]])) # old A
    } else {
      C_tilde[[j]] <- t(twodG(A = Y@B[[j]] %*% Y@C[[j]], B = Y@B[[j]], grid = Y@grid[[j]]))
      G_1[[j]] <- t(twodG(A = Y@B[[j]], B = Y@B[[j]], grid = Y@grid[[j]]))
    }
    shifter[1L, j + 1L] <- shifter[2L, j] + 1L
    shifter[2L, j + 1L] <- shifter[2L, j] + L*ncol(Y@B[[j]])
    Y_d[1,j] = ncol(Y@B[[j]])
  }
  d <- cbind(0,Y_d*L)
  # Find the proper inner product matrices for j_k variables
  d_tilde <- sum(d) / L
  K <- N - L + 1L
  # Calculating Variance/Covariance Structure
  S0 <- SSM(K, L, d_tilde, p, C_tilde, shifter)
  # Calculating Gram Matrix
  H <- CalculateInverse(Gramm(K, L, p, d_tilde, G_1, shifter, d))
  # Calculating Eigen Triples
  Q <- eigs(AtimesB(H, S0), ntriples)
  # Returning results
  Q$values <- Re(Q$values)
  Q$vectors <- Re(Q$vectors)
  coefs0 <- Q$vectors
  p_c <- list()
  values <- Q$values[1L:ntriples]
  out <- list()
  for (i in 1L:(ntriples)) {
    my_pcs <- list(NA)
    for (j in 1L:p) {
      my_pcs[[j]] <- Y@B[[j]] %*% Cofmat((d[j + 1L] / L), L, coefs0[(shifter[1L, (j + 1L)]:shifter[2L, (j + 1L)]), i])
    }
    out[[i]] <- my_pcs
  }
  out$values <- values
  out$L <- L
  out$N <- N
  out$Y <- Y
  out$RVectrs <- mV(out, ntriples)
  return(out)
}
