#--------------------------------------------------------------
#' Reconstruction of Functional TS
#'
#' is a function for reconstruction stage (including Grouping and
#' Hankelization steps) The output is a list of functional time series corresponds to each group.
#' 'U' in the input is a fssa object. 'group' is a list.
#' @param U funtional object
#' @return matrix contating Cx
#' @export

freconstruct <- function(U, group = as.list(1L:10L)) {
  N <- U$N
  Y <- U$Y
  d <- nrow(U[[1L]]$coefs)
  L <- U$L
  K <- N - L + 1L
  basis <- U[[1]]$basis
  m <- length(group)
  basis <- U[[1L]]$basis
  out <- list()
  for (i in 1L:m) {
    Cx <- matrix(NA, nr = d, nc = N)
    g <- group[[i]]
    S <- 0L
    for (j in 1L:length(g)) S <- S +
      fproj(U, g[j], d, K, L, Y)
    S <- fH(S, d)
    Cx[, 1L:L] <- S[, 1L, ]
    Cx[, L:N] <- S[, ,L]
    out[[i]] <- fd(Cx, basis)
  }
  out$values <- sqrt(U$values)
  return(out)
}
