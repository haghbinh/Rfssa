# Embedding and decomposition stages of multivariate functional singular spectrum analysis.
mfssa <- function(Y, L, ntriples) {
  # get c plus plus code
  p <- length(Y@C)
  Y_d <- matrix(data = 0, nrow = 1, ncol = p)
  for (j in 1:p) Y_d[1, j] <- ncol(Y@B[[j]])
  d <- L * matrix(c(0, Y_d), nrow = 1L, ncol = (p + 1L))
  N <- ncol(Y@C[[1]])
  C_tilde <- list()
  G_1 <- list()
  # get inner product matrices
  for (i in 1:p) {
    if (ncol(Y@grid[[i]]) == 1) {
      C_tilde[[i]] <- t(onedG(A = Y@B[[i]] %*% Y@C[[i]], B = Y@B[[i]], grid = Y@grid[[i]])) # old B
      G_1[[i]] <- t(onedG(A = Y@B[[i]], B = Y@B[[i]], grid = Y@grid[[i]])) # old A
    } else {
      C_tilde[[i]] <- t(twodG(A = Y@B[[i]] %*% Y@C[[i]], B = Y@B[[i]], grid = Y@grid[[i]]))
      G_1[[i]] <- t(twodG(A = Y@B[[i]], B = Y@B[[i]], grid = Y@grid[[i]]))
    }
  }
  # Find the proper inner product matrices for j_k variables
  d_tilde <- sum(d) / L
  K <- N - L + 1L
  shifter <- matrix(nrow = 2, ncol = (p + 1L), data = 0L)
  shifter[, 2L] <- c(1L, d[2L])
  if (p > 1L) {
    for (i in 2L:p) {
      shifter[1L, i + 1L] <- shifter[2L, i] + 1L
      shifter[2L, i + 1L] <- shifter[2L, i] + d[i + 1L]
    }
  }
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
