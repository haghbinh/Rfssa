#' Functional SSA
#'
#' fssa is a function for decomposition stage (including embeding
#'  and  SVD step) of a functional time series.
#' @return list The outputs of following function is a list which
#' includes functional eigen vectors in the form of L-variate functional object.
#' @param Y A functional time series.
#' @param L Windows length
#' @importFrom fda fd
#' @importFrom fda inprod eval.fd smooth.basis
#' @useDynLib Rfssa
#' @export
fssa <- function(Y, L = floor(dim(Y$coefs)[2L] / 2L)) {
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
    class(out) <- "fssa"
    return(out)
}

