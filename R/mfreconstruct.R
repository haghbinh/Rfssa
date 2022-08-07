# Reconstruction stage (including Hankelization) of multivariate functional singular spectrum analysis.
mfreconstruct <- function(U, groups = as.list(1L:10L)) {
  N <- U$N
  Y <- U$Y
  p <- length(U$Y@C)
  L <- U$L
  K <- N - L + 1
  basis <- U$Y@B
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
      if (ncol(Y@grid[[j]]) == 2) {
        x <- unique(Y@grid[[j]][, 1])
        y <- unique(Y@grid[[j]][, 2])
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
        new_grid[[j]] <- Y@grid[[j]]
      }
    }



    # output the reconstructions
    fts_out <- Rfssa::fts(X = recon_group, B = basis, grid = new_grid,time = colnames(U$Y@C[[1]]))
    fts_out@basis_type <- Y@basis_type
    recon_out[[i]] <- fts_out
  }
  recon_out$values <- U$values
  return(recon_out)
}
