# FSSA Recurrent and Vector Forecasting of Multivariate FTS

mfforecast <- function(U, groups = list(c(1)), h = 1, method = "recurrent", tol = 10^-3) {
  out <- list()
  if (method == "vector") {
    warning("MFSSA vector forecasting displays some instabilities for certain datasets and certain chosen forecast horizons. The instability phenomenon is still under investigation. Use the MFSSA vector forecasting algorithm with caution or set method=\"recurrent\" to use the stable MFSSA recurrent forecasting technique.")
  }
  for (a in 1:length(groups)) {
    g <- groups[[a]]
    # Define prediction space
    basis <- list()
    p <- length(U$Y@C)
    d <- matrix(data = 0, nrow = 1, ncol = p)
    N <- U$N
    L <- U$L
    K <- N - L + 1
    shifter <- matrix(data = 0, nrow = 2, ncol = p)
    shifter[1, 1] <- 1
    shifter[2, 1] <- ncol(U$Y@B[[1]])
    k <- length(g)
    D <- matrix(data = NA, nrow = 1, ncol = k)
    for (j in 1:p) {
      basis[[j]] <- U$Y@B[[j]]
      d[j] <- ncol(U$Y@B[[j]])
      if (j > 1) {
        shifter[1, j] <- shifter[2, (j - 1)] + 1
        shifter[2, j] <- shifter[2, (j - 1)] + ncol(U$Y@B[[j]])
      }
      D_j <- matrix(data = NA, nrow = d[j], ncol = k)
      for (n in 1:k) D_j[, n] <- solve(t(basis[[j]]) %*% basis[[j]]) %*% t(basis[[j]]) %*% U[[g[n]]][[j]][, L]
      D <- rbind(D, D_j)
    }
    D <- D[(2:nrow(D)), ]

    # Define Truncated Neumann Series
    NU_1 <- D %*% t(D)
    NU <- list()
    NU[[1]] <- diag(1, sum(d))
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
      # MFSSA R-forecasting
      out_g <- list()
      Q <- freconstruct(U, groups = list(g))
      fssa_fore <- matrix(data = 0, nrow = sum(d), ncol = h)
      for (m in 1:h) {
        for (j in 1:(L - 1)) {
          my_obs <- matrix(data = 0, nrow = sum(d), ncol = 1)
          E_j <- matrix(data = NA, nrow = sum(d), ncol = k)
          for (q in 1:p) {
            for (n in 1:k) {
              E_j[(shifter[1, q]:shifter[2, q]), n] <- solve(t(basis[[q]]) %*% basis[[q]]) %*% t(basis[[q]]) %*% U[[g[n]]][[q]][, j]
            }
            my_obs[(shifter[1, q]:shifter[2, q]), 1] <- Q[[1]]@C[[q]][, (N + j - L + m)]
          }
          A_j <- Neu %*% D %*% t(E_j)
          fssa_fore[, m] <- fssa_fore[, m] + A_j %*% my_obs
        }

        for (q in 1:p) {
          Q[[1]]@C[[q]] <- cbind(Q[[1]]@C[[q]], fssa_fore[(shifter[1, q]:shifter[2, q]), m])
        }
      }
      for (q in 1:p) {
        out_g[[q]] <- basis[[q]] %*% fssa_fore[(shifter[1, q]:shifter[2, q]), ]
      }
      out[[a]] <- Rfssa::fts(out_g, basis, U$Y@grid)
    } else if (method == "vector") {
      # MFSSA V-forecasting
      out_g <- list()
      Y <- matrix(data = NA, nrow = ((L - 1) * sum(d)), ncol = k)
      Lshift <- shifter
      Lshift[2, 1] <- (L - 1) * d[1]
      for (j in 2:p) Lshift[1, j] <- (Lshift[2, (j - 1)]) + 1
      Lshift[2, j] <- Lshift[2, (j - 1)] + (L - 1) * d[j]
      for (j in 1:p) {
        for (n in 1:k) {
          Y[(Lshift[1, j]:Lshift[2, j]), n] <- matrix(data = solve(t(basis[[j]]) %*% basis[[j]]) %*% t(basis[[j]]) %*% U[[g[n]]][[j]][, 1:(L - 1)], nrow = ((L - 1) * d[j]), ncol = 1)
        }
      }
      P <- Y %*% t(Y) + Y %*% t(D) %*% Neu %*% D %*% t(Y)
      S <- list()
      for (j in 1:p) {
        S[[j]] <- array(data = 0, dim = c(d[j], (K + h), L))
      }
      for (k in 1L:length(g)) {
        projection <- mfproj(U, g[k])
        for (j in 1:p) {
          S[[j]][, 1:K, ] <- S[[j]][, 1:K, ] + projection[[j]]
        }
      }
      for (m in 1:h) {
        obs <- matrix(data = NA, nrow = ((L - 1) * sum(d)), ncol = 1)
        for (j in 1:p) {
          obs[(Lshift[1, j]:Lshift[2, j]), 1] <- matrix(data = S[[j]][, (K + (m - 1)), 2:L], nrow = ((L - 1) * d[j]), ncol = 1)
        }
        pr <- P %*% obs
        pr_mat <- matrix(data = NA, nrow = sum(d), ncol = L - 1)
        for (j in 1:p) {
          pr_mat[(shifter[1, j]:shifter[2, j]), ] <- pr[(Lshift[1, j]:Lshift[2, j]), 1]
        }
        pr_1 <- matrix(data = 0, nrow = sum(d), ncol = 1)
        for (j in 1:(L - 1)) {
          my_obs <- matrix(data = 0, nrow = sum(d), ncol = 1)
          E_j <- matrix(data = NA, nrow = sum(d), ncol = k)
          for (q in 1:p) {
            for (n in 1:k) {
              E_j[(shifter[1, q]:shifter[2, q]), n] <- solve(t(basis[[q]]) %*% basis[[q]]) %*% t(basis[[q]]) %*% U[[g[n]]][[q]][, j]
            }
            my_obs <- pr_mat[, j]
          }
          A_j <- Neu %*% D %*% t(E_j)
          pr_1 <- pr_1 + A_j %*% my_obs
        }


        for (q in 1:p) {
          S[[q]][, (K + m), ] <- matrix(data = cbind(pr_mat, pr_1)[(shifter[1, q]:shifter[2, q]), ], nrow = d[q], ncol = L)
        }
      }

      out_g <- list()
      for (q in 1:p) {
        S[[q]] <- fH(S[[q]], d[q])
        out_g[[q]] <- basis[[q]] %*% S[[q]][, (K + 1):(K + h), L]
      }


      out[[a]] <- Rfssa::fts(out_g, basis, U$Y@grid)
    }
  }

  return(out)
}
