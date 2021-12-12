# Embedding and decomposition stages of univariate functional singular spectrum analysis
ufssa <- function(Y, L, ntriples) {
  N <- ncol(Y@C[[1]])
  basis <- Y@B[[1]]
  d <- ncol(Y@B[[1]])
  grid <- Y@grid[[1]]
  K <- N - L + 1L
  if (ncol(Y@grid[[1]]) == 1) {
    C_tilde <- t(onedG(A = basis %*% Y@C[[1]], B = basis, grid = Y@grid[[1]]))
    G <- onedG(A = basis, B = basis, grid = Y@grid[[1]])
  } else {
    C_tilde <- t(twodG(A = basis %*% Y@C[[1]], B = basis, grid = Y@grid[[1]]))
    G <- twodG(A = basis, B = basis, grid = Y@grid[[1]])
  }
  # Calculating Variance/Covariance Structure
  S0 <- SS(K, L, C_tilde, d)
  # Calculating Gram Matrix
  H <- CalculateInverse(Gram(K, L, G, d))
  # Calculating Eigen Triples
  Q <- eigs(AtimesB(H, S0), ntriples)
  # Returning results
  Q$values <- Re(Q$values)
  Q$vectors <- Re(Q$vectors)
  out <- list(NA)
  for (i in 1L:ntriples) out[[i]] <- Y@B[[1]] %*% Cofmat(d, L, Q$vectors[, i])
  out$values <- Q$values[1L:ntriples]
  out$L <- L
  out$N <- N
  out$Y <- Y
  out$RVectrs <- uV(out, ntriples)
  return(out)
}
