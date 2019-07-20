# Univariate FSSA
ufssa <- function(Y, L) {
  N <- dim(Y$coefs)[2]
  basis <- Y$basis
  d <- basis$nbasis
  K <- N - L + 1L
  B <- inprod(Y, basis)
  A <- inprod(basis, basis)
  S0 <- SS(K, L, B, d)
  H <- solve(Gram(K, L, A, d))
  Q <- eigen(H %*% S0)
  Q$vectors <- Re(Q$vectors)
  out <- list(NA)
  d1 <- sum(Re(Q$values) > 0.001)
  for (i in 1L:(d1)) out[[i]] <- fd(Cofmat(d, L, Q$vectors[, i]), basis)
  out$values <- Re(Q$values[1L:d1])
  out$L <- L
  out$N <- N
  out$Y <- Y
}

