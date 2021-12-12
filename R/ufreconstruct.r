# Reconstruction stage (including Hankelization) of univariate functional singular spectrum analysis

ufreconstruct <- function(U, groups = as.list(1L:10L)) {
  N <- ncol(U$Y@C[[1]])
  Y <- U$Y
  d <- ncol(U$Y@B[[1]])
  L <- U$L
  K <- N - L + 1L
  basis <- Y@B[[1]]
  m <- length(groups)
  out <- list()
  for (i in 1L:m) {
    Cx <- matrix(NA, nrow = d, ncol = N)
    g <- groups[[i]]
    S <- 0L
    for (j in 1L:length(g)) S <- S + ufproj(U, g[j], d)
    S <- fH(S, d)
    Cx[, 1L:L] <- S[, 1L, ]
    Cx[, L:N] <- S[, , L]
    recon_out <- basis %*% Cx
    if (ncol(Y@grid[[1]]) == 2) {
      x <- unique(Y@grid[[1]][, 1])
      y <- unique(Y@grid[[1]][, 2])
      recon_two_d <- array(data = NA, dim = c(length(x), length(y), N))
      for (n in 1:N) {
        count <- 1
        for (i_1 in 1:length(x)) {
          for (i_2 in 1:length(y)) {
            recon_two_d[i_1, i_2, n] <- recon_out[count, n]
            count <- count + 1
          }
        }
      }

      recon_out <- recon_two_d
      new_grid <- list(x, y)
    } else {
      new_grid <- Y@grid[[1]]
    }
    out[[i]] <- Rfssa::fts(list(recon_out), list(Y@B[[1]]), list(new_grid))
  }
  out$values <- sqrt(U$values)
  return(out)
}
