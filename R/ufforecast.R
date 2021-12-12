ufforecast <- function(U, groups = list(c(1)), h = 1, method = "recurrent", tol = 10^-3) {
  out <- list()
  for (a in 1:length(groups)) {
    g <- groups[[a]]
    # Define prediction space
    basis <- U$Y@B[[1]]
    N <- U$N
    d <- ncol(basis)
    L <- U$L
    K <- N - L + 1
    k <- length(g)
    G <- t(basis) %*% basis
    D <- matrix(data = NA, nrow = d, ncol = k)
    for (n in 1:k) D[, n] <- solve(G) %*% t(basis) %*% U[[g[n]]][, L]

    # Define Truncated Neumann Series
    NU_1 <- D %*% t(D)
    NU <- list()
    NU[[1]] <- diag(1, d)
    norm <- 1
    l <- 0
    while (norm > tol) {
      l <- l + 1
      norm <- norm(NU[[l]])
      NU[[(l + 1)]] <- NU_1 %*% NU[[l]]
    }
    Neu <- NU[[1]]
    for (n in 2:l) Neu <- Neu + NU[[n]]


    if (method == "recurrent") {
      # FSSA R-forecasting

      # Reconstruct signal
      Q <- Rfssa::freconstruct(U, groups = list(g))
      fssa_fore <- matrix(data = 0, nrow = d, ncol = h)
      for (m in 1:h) {
        for (j in 1:(L - 1)) {
          E_j <- matrix(data = NA, nrow = d, ncol = k)
          for (n in 1:k) {
            E_j[, n] <- solve(G) %*% t(basis) %*% U[[g[n]]][, j]
          }
          A_j <- Neu %*% D %*% t(E_j)
          fssa_fore[, m] <- fssa_fore[, m] + A_j %*% as.matrix(Q[[1]]@C[[1]][, (N + j - L + m)])
        }
        Q[[1]]@C[[1]] <- cbind(Q[[1]]@C[[1]], fssa_fore[, m])
      }
      out[[a]] <- Rfssa::fts(list(basis %*% fssa_fore), list(basis), list(Q[[1]]@grid[[1]]))
    } else if (method == "vector") {
      # FSSA V-forecasting
      F_matrix <- matrix(data = NA, nrow = ((L - 1) * d), ncol = k)
      for (n in 1:k) {
        F_matrix[, n] <- matrix(data = solve(G) %*% t(basis) %*% U[[g[n]]][, 1:(L - 1)], nrow = ((L - 1) * d), ncol = 1)
      }
      P <- F_matrix %*% t(F_matrix) + F_matrix %*% t(D) %*% Neu %*% D %*% t(F_matrix)
      S <- array(data = 0, dim = c(d, (K + h), L))
      for (j in 1L:length(g)) {
        S[, (1:K), ] <- S[, (1:K), ] + ufproj(U, g[j], d)
      }
      for (m in 1:h) {
        obs <- matrix(data = S[, (K + (m - 1)), 2:L], nrow = ((L - 1) * d), ncol = 1)
        pr <- P %*% obs
        pr <- matrix(data = pr, nrow = d, ncol = (L - 1))
        pr_1 <- matrix(data = 0, nrow = d, ncol = 1)
        for (j in 1:(L - 1)) {
          E_j <- matrix(data = NA, nrow = d, ncol = k)
          for (n in 1:k) {
            E_j[, n] <- solve(G) %*% t(basis) %*% U[[g[n]]][, j]
          }
          A_j <- Neu %*% D %*% t(E_j)

          pr_1[, 1] <- pr_1[, 1] + A_j %*% S[, (K + m - 1), (j + 1)]
        }
        pr <- cbind(pr, pr_1)
        S[, (K + m), ] <- pr
      }
      S <- fH(S, d)
      predictions <- S[, (K + 1):(K + h), L]

      out[[a]] <- Rfssa::fts(list(basis %*% predictions), list(basis), list(U$Y@grid[[1]]))
    }
  }
  return(out)
}
