# =============================================================================
# Rank selection for PrSSA by composite-likelihood information criteria.
# =============================================================================

#' Rank Selection for Probabilistic SSA
#'
#' Fits \code{\link{prssa}} for \eqn{r = 1, \dots, r_{max}} and compares the fits by
#' composite-likelihood information criteria that account for the overlap of the windows.
#'
#' @details
#' With per-window scores \eqn{s_j} and the Fisher information \eqn{H} of one window, the
#' effective number of parameters is \eqn{p^{eff}_r = \mathrm{tr}(\hat J \hat H^{+})}, where
#' \eqn{\hat J} is a Bartlett HAC estimate of the long-run variance of the scores with
#' bandwidth \eqn{M = L - 1 + \lfloor 4(K/100)^{2/9}\rfloor} (overlapping windows always have
#' correlated scores, so \eqn{M \ge L - 1}). The criteria are
#' \itemize{
#'   \item \code{clbic} (default): \eqn{-2\ell_C + p^{eff}_r \log K};
#'   \item \code{clbic_nominal}: \eqn{-2\ell_C + p_r \log K} with the nominal count
#'   \eqn{p_r = Dr - r(r-1)/2 + } (number of noise parameters); it needs only the fits;
#'   \item \code{claic}: \eqn{-2\ell_C + 2 p^{eff}_r} (tends to overfit).
#' }
#' The trace is computed exactly without forming the full information matrix: within each
#' isotropic noise block, the directions orthogonal to the block's loadings decouple, so only a
#' system of size \eqn{S r + B} (\eqn{S \le B r}) is inverted. For heteroscedastic noise
#' (\code{noise = "hetero"}) the sandwich penalty is not available; the nominal penalty is
#' then used for \code{clbic} and \code{claic}, with a warning.
#'
#' @param Y,L Data and window length, as in \code{\link{prssa}}.
#' @param rmax Largest rank considered.
#' @param criterion Criterion used to select the rank: \code{"clbic"} (default),
#'   \code{"clbic_nominal"} or \code{"claic"}.
#' @param ... Further arguments passed to \code{\link{prssa}} (\code{noise}, \code{mean},
#'   \code{alpha}, ...).
#'
#' @return An object of class \code{"prrank"}: a list with \code{table} (a data frame with
#'   columns \code{r}, \code{loglik}, \code{p_nominal}, \code{p_eff}, \code{clbic},
#'   \code{clbic_nominal}, \code{claic}), \code{selected} (the selected rank) and
#'   \code{criterion}.
#'
#' @references
#' Varin, C., Reid, N., & Firth, D. (2011). An overview of composite likelihood methods.
#' \emph{Statistica Sinica}, \strong{21}(1), 5-42.
#'
#' @seealso \code{\link{prssa}}
#'
#' @examples
#' set.seed(1)
#' tt <- 1:200
#' y <- 2 * sin(2 * pi * tt / 20) + cos(2 * pi * tt / 7) + rnorm(200, sd = 0.7)
#' rk <- prrank(y, L = 24, rmax = 6)
#' rk
#'
#' @export
prrank <- function(Y, L, rmax, criterion = c("clbic", "clbic_nominal", "claic"), ...) {
  criterion <- match.arg(criterion)
  if (missing(rmax) || length(rmax) != 1 || rmax < 1 || rmax != round(rmax)) {
    stop("`rmax` must be a positive integer.")
  }
  tab <- NULL
  warned <- FALSE
  for (r in seq_len(rmax)) {
    fit <- prssa(Y, L, r, ...)
    pe <- if (fit$noise == "hetero") NA_real_ else .pr_peff(fit)
    if (is.na(pe) && !warned && criterion != "clbic_nominal") {
      warning("The sandwich penalty is not available for noise = \"hetero\"; the nominal parameter count is used.")
      warned <- TRUE
    }
    pen <- if (is.na(pe)) fit$df else pe
    tab <- rbind(tab, data.frame(
      r = fit$r, loglik = fit$loglik, p_nominal = fit$df, p_eff = pe,
      clbic = -2 * fit$loglik + pen * log(fit$K),
      clbic_nominal = -2 * fit$loglik + fit$df * log(fit$K),
      claic = -2 * fit$loglik + 2 * pen
    ))
  }
  sel <- tab$r[which.min(tab[[criterion]])]
  structure(list(table = tab, selected = sel, criterion = criterion), class = "prrank")
}

# Bartlett HAC weights w_1, ..., w_M with M = L - 1 + floor(4 (K/100)^(2/9))
.pr_hac_weights <- function(L, K) {
  M <- min(L - 1 + floor(4 * (K / 100)^(2 / 9)), K - 1)
  if (M < 1) return(numeric(0))
  1 - seq_len(M) / (M + 1)
}

# HAC trace  sum_m c_m (1/K) sum_j a_j b_{j+m}  for per-window K x p matrices
.pr_hac_trace <- function(A, B, w) {
  K <- nrow(A)
  tr <- sum(A * B) / K
  for (m in seq_along(w)) {
    tr <- tr + w[m] * (sum(A[1:(K - m), , drop = FALSE] * B[(1 + m):K, , drop = FALSE]) +
                         sum(B[1:(K - m), , drop = FALSE] * A[(1 + m):K, , drop = FALSE])) / K
  }
  tr
}

# Moore-Penrose inverse of a symmetric matrix
.pr_pinv <- function(A, tol = 1e-9) {
  e <- eigen((A + t(A)) / 2, symmetric = TRUE)
  k <- e$values > tol * max(e$values)
  e$vectors[, k, drop = FALSE] %*% (t(e$vectors[, k, drop = FALSE]) / e$values[k])
}

# effective number of parameters tr(J H^+) for block-isotropic noise (exact,
# reduced: noise directions orthogonal to the block loadings decouple)
.pr_peff <- function(fit) {
  X <- .pr_fit_trajectory(fit, smoothed = TRUE)
  W <- if (fit$alpha > 0) {
    sqrt(fit$psi) * sweep(fit$U, 2, sqrt(fit$lambda - 1), "*")
  } else {
    fit$W
  }
  psi <- fit$psi
  K <- ncol(X)
  D <- nrow(X)
  r <- ncol(W)
  blk <- if (fit$nz$iso) rep(1L, D) else fit$nz$unit
  B <- max(blk)
  rows_b <- split(seq_len(D), blk)
  sig2 <- vapply(rows_b, function(rr) psi[rr[1]], numeric(1))
  Wp <- W / psi
  A0 <- crossprod(W, Wp)
  Minv <- solve(diag(r) + A0)
  CiE <- X / psi - Wp %*% (Minv %*% crossprod(Wp, X))
  P <- Wp - Wp %*% (Minv %*% A0)
  Q <- A0 - A0 %*% Minv %*% A0
  Qinv <- solve(Q)
  G <- crossprod(W, CiE)
  wts <- .pr_hac_weights(fit$L, K)
  # noise part: sum_b eps_jb' eps_kb / sigma_b^2 times g_j' Q^{-1} g_k
  Tb <- lapply(rows_b, function(rr) .pr_orth(W[rr, , drop = FALSE]))
  Eps2 <- matrix(0, 0, K)
  a_rows <- list()
  for (b in seq_len(B)) {
    rr <- rows_b[[b]]
    Xb <- X[rr, , drop = FALSE]
    Eb <- Xb - Tb[[b]] %*% crossprod(Tb[[b]], Xb)
    a_rows[[b]] <- Eb / sqrt(sig2[b])
  }
  Ea <- do.call(rbind, a_rows)
  Gq <- .pr_msqrt(Qinv) %*% G
  noise_tr <- {
    pair <- function(m) sum(colSums(Ea[, 1:(K - m), drop = FALSE] * Ea[, (1 + m):K, drop = FALSE]) *
                              colSums(Gq[, 1:(K - m), drop = FALSE] * Gq[, (1 + m):K, drop = FALSE])) / K
    tr <- sum(colSums(Ea^2) * colSums(Gq^2)) / K
    for (m in seq_along(wts)) tr <- tr + 2 * wts[m] * pair(m)
    tr
  }
  # signal part: parameters T'W (S x r) and the B noise variances
  Tm <- matrix(0, D, 0)
  for (b in seq_len(B)) {
    Tb_full <- matrix(0, D, ncol(Tb[[b]]))
    Tb_full[rows_b[[b]], ] <- Tb[[b]]
    Tm <- cbind(Tm, Tb_full)
  }
  S <- ncol(Tm)
  TCiE <- crossprod(Tm, CiE)
  PS <- crossprod(Tm, P)
  SW <- do.call(cbind, lapply(seq_len(r), function(a) t(TCiE) * G[a, ])) - rep(as.numeric(PS), each = K)
  trCi <- vapply(rows_b, function(rr) length(rr) / psi[rr[1]] - sum((Wp[rr, , drop = FALSE] %*% Minv) * Wp[rr, , drop = FALSE]), numeric(1))
  Ss <- vapply(seq_len(B), function(b) 0.5 * (colSums(CiE[rows_b[[b]], , drop = FALSE]^2) - trCi[b]), numeric(K))
  if (K == 1) Ss <- matrix(Ss, 1)
  Smat <- cbind(SW, Ss)
  TCiT <- diag(rep(1 / sig2, vapply(Tb, ncol, integer(1))), S) - crossprod(Tm, Wp) %*% Minv %*% crossprod(Wp, Tm)
  iv <- rep(seq_len(S), r)
  bv <- rep(seq_len(r), each = S)
  F1 <- PS[iv, bv, drop = FALSE]
  HWW <- F1 * t(F1) + kronecker(Q, TCiT)
  HWs <- vapply(seq_len(B), function(b) {
    EbP <- matrix(0, D, r)
    EbP[rows_b[[b]], ] <- P[rows_b[[b]], , drop = FALSE]
    CiEbP <- EbP / psi - Wp %*% (Minv %*% crossprod(Wp, EbP))
    as.numeric(crossprod(Tm, CiEbP))
  }, numeric(S * r))
  if (S * r == 1) HWs <- matrix(HWs, 1)
  Hss <- matrix(0, B, B)
  for (b in seq_len(B)) {
    for (cc in seq_len(B)) {
      Wb <- Wp[rows_b[[b]], , drop = FALSE]
      Wc <- Wp[rows_b[[cc]], , drop = FALSE]
      val <- sum(diag(Minv %*% crossprod(Wb) %*% Minv %*% crossprod(Wc)))
      if (b == cc) {
        val <- val + length(rows_b[[b]]) / sig2[b]^2 -
          2 * sum((Wb %*% Minv) * Wb) / sig2[b]
      }
      Hss[b, cc] <- 0.5 * val
    }
  }
  H <- rbind(cbind(HWW, HWs), cbind(t(HWs), Hss))
  Hp <- .pr_pinv(H)
  signal_tr <- .pr_hac_trace(Smat, Smat %*% Hp, wts)
  noise_tr + signal_tr
}
