# =============================================================================
# PrSSA reconstruction: posterior-mean (shrunken SSA) or projection records,
# followed by diagonal averaging.
# =============================================================================

#' Reconstruction for Probabilistic SSA
#'
#' Reconstructs the signal of a \code{\link{prssa}} fit by grouping its components. Each
#' lagged record is replaced by the posterior mean of its signal part (\code{"posterior"}) or
#' by its projection on the signal subspace (\code{"projection"}), and the overlapping
#' estimates are combined by diagonal averaging.
#'
#' @details
#' With noise-whitened singular triples \eqn{(\tilde s_k, \tilde u_k, \tilde v_k)} of the
#' trajectory matrix, component \eqn{k} contributes
#' \eqn{\omega_k \Psi^{1/2}\tilde u_k \tilde s_k \tilde v_k^\top}, where
#' \eqn{\omega_k = 1 - 1/\tilde\lambda_k} (posterior mean) or \eqn{\omega_k = 1} (projection).
#' For isotropic noise \eqn{\omega_k = 1 - \hat\sigma^2/\kappa_k} and the projection
#' reconstruction is exactly the SSA (FSSA, VMSSA) reconstruction of the same group; the
#' posterior reconstruction is SSA with component-specific shrinkage. For block noise the
#' projection reconstruction is weighted MSSA with weights \eqn{\hat\sigma_b^{-2}}. With
#' smoothing (\code{alpha > 0}) the half-smoothing operator is applied to the result
#' (regularized HSSA in the noise-free limit).
#'
#' The reconstructions do not contain the mean; it is returned separately in \code{$mean}, so
#' that \code{mean + } (reconstruction of group \code{1:r}) estimates the signal.
#'
#' @param fit An object of class \code{"prssa"}.
#' @param groups A list of integer vectors, each a group of component indices in
#'   \code{1:fit$r}; a single vector is treated as one group.
#' @param type \code{"posterior"} (default) or \code{"projection"}.
#'
#' @return A list with one element per group (named \code{F1}, \code{F2}, ... unless
#'   \code{groups} is named), each of the same class as the input series (numeric vector or
#'   \code{ts}, matrix, list, or \code{\link{funts}}), followed by \code{mean} (the fitted mean
#'   in the same class) and \code{omega} (the shrinkage factors used).
#'
#' @seealso \code{\link{prssa}}, \code{\link{freconstruct}}
#'
#' @examples
#' set.seed(1)
#' tt <- 1:200
#' y <- 2 * sin(2 * pi * tt / 20) + cos(2 * pi * tt / 7) + rnorm(200, sd = 0.7)
#' fit <- prssa(y, L = 24, r = 4)
#' rec <- prreconstruct(fit, groups = list(1:2, 3:4))
#' plot(tt, y, col = "grey", pch = 16, cex = 0.5)
#' lines(tt, rec$mean + rec$F1 + rec$F2, lwd = 2)
#'
#' @export
prreconstruct <- function(fit, groups = list(1:fit$r), type = c("posterior", "projection")) {
  if (!inherits(fit, "prssa")) stop("`fit` must be an object of class \"prssa\".")
  type <- match.arg(type)
  if (is.numeric(groups) && !is.list(groups)) groups <- list(groups)
  allidx <- unlist(groups)
  if (length(allidx) == 0 || any(!allidx %in% seq_len(fit$r))) {
    stop(sprintf("Group indices must lie in 1:%d (the retained components).", fit$r))
  }
  dat <- fit$dat
  L <- fit$L
  dims <- dat$dims
  omega <- if (type == "posterior") fit$omega else rep(1, fit$r)
  sq <- sqrt(fit$psi)
  Slist <- lapply(dat$blocks, function(b) b$S)
  out <- vector("list", length(groups))
  for (i in seq_along(groups)) {
    g <- groups[[i]]
    Tm <- (sq * fit$U[, g, drop = FALSE]) %*% ((omega[g] * fit$s[g]) * t(fit$V[, g, drop = FALSE]))
    if (fit$alpha > 0) Tm <- .pr_apply_blockdiag(Tm, Slist, L, dims)
    out[[i]] <- .pr_as_input(dat, .pr_hankelize(Tm, L, dims), seq_len(fit$N))
  }
  names(out) <- if (!is.null(names(groups))) names(groups) else paste0("F", seq_along(groups))
  mser <- lapply(fit$mean, function(m) matrix(m, length(m), fit$N))
  out$mean <- .pr_as_input(dat, mser, seq_len(fit$N))
  out$omega <- omega
  out
}
