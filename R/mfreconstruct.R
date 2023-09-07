# Reconstruction stage (including Hankelization) of multivariate functional singular spectrum analysis.
mfreconstruct <- function(U, groups = as.list(1L:10L)) {
  Y <- U$Y
  N <- Y$N
  basisobj <- Y$basis
  time_st <- Y$time[1]
  time_en <- Y$time[N]
  L <- U$L
  K <- N - L + 1
  basis <- Y$B_mat
  m <- length(groups)
  recon_out <- list()
  new_grid <- list()

  # Loop over groups
  for (i in 1L:m) {
    recon_group <- list()
    C <- list()
    S <- list()
    g <- groups[[i]]
    # build reconstructions
    for (j in 1:p) {
      d <- ncol(basis[[j]])
      C[[j]] <- matrix(NA, nrow = d, ncol = N)
      S[[j]] <- 0L
      for(k in 1L:length(g)){
        S[[j]] <- S[[j]] + mfproj(U, g[k])[[j]]
      }
      S[[j]] <- fH(S[[j]], d)
      C_jx <- C[[j]]
      S_jx <- S[[j]]
      C_jx[, 1L:L] <- S_jx[, 1L, ]
      C_jx[, L:N] <- S_jx[, , L]
      recon_group[[j]] <- basis[[j]] %*% C_jx
      if (ncol(Y$argval[[j]]) == 2) {
        x <- unique(Y$argval[[j]][, 1])
        y <- unique(Y$argval[[j]][, 2])
        recon_two_d <- array(data = NA, dim = c(length(x), length(y), N))
        for (n in 1:N) {
          count <- 1
          for (i_1 in 1:length(x)) {
            for (i_2 in 1:length(y)) {
              recon_two_d[i_1, i_2, n] <- recon_group[[j]][count, n]
              count <- count + 1
            }
          }
        }

        recon_group[[j]] <- recon_two_d
        new_grid[[j]] <- list(x, y)
      } else {
        new_grid[[j]] <- Y$argval[[j]]
      }
    }

    # output the reconstructions
    funts_out <- funts(X = recon_group, basisobj = basisobj, argval = new_grid, method = "coefs", start = time_st, end = time_en)
    recon_out[[i]] <- fts_out
  }
  recon_out$values <- U$values
  return(recon_out)
}
