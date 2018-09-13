#' Wcorrelation matrix
#'
#' This function evaluate the Wcorrolation plot for fssa
#' @param U in the input is a fssa object.
#' @param d is the number of elementary components.
#' @export
fwcor <- function(U, d) {
  Q <- freconstruct(U, group = as.list(1:d))
  N <- U$N
  L <- ncol(U[[1]]$coefs)
  K <- N - L + 1L
  w <- 1L:N
  basis <- U[[1]]$basis
  G <- inprod(basis, basis)
  L1 <- min(L, K)
  K1 <- max(K, L)
  w[L1:K1] <- L1
  w[(K1 + 1L):N] <- N + 1L - ( (K1 +1L):N)
  out <- matrix(1L, nrow = d, ncol = d)
  for (i in 1L:(d - 1)) {
    for (j in (i + 1L):d) {
      out[i, j] <- winprod(Q[[i]], Q[[j]], w, G)/sqrt(winprod(Q[[i]],Q[[i]], w, G) * winprod(Q[[j]],Q[[j]], w, G))
    }
  }
  for (i in 2:d) for (j in 1:(i - 1)) out[i, j] <- out[j, i]
  return(out)
}


