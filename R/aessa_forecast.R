# ------------------------------------------------------------------------------
# Forecasting for AESSA / VAESSA: decoder inversion (window conditioning),
# recursion, and latent dynamics.
# ------------------------------------------------------------------------------

# Levenberg-Marquardt solution of
#   min_h  sum_o w_o (x_o - g_o(h))^2 + lambda ||h||^2
# over the observed record coordinates O, from every row of `starts`.
.ae_invert <- function(par, spec, xo, O, w, lambda, starts, tg = NULL, maxit = 100L, xtol = 1e-10) {
  r <- spec$r
  sw <- sqrt(w)
  eval_at <- function(z) {
    jj <- .ae_jac_R(par, z, spec)
    res <- xo - jj$g[O]
    list(jj = jj, res = res, f = sum(w * res^2) + lambda * sum(z^2))
  }
  sols <- vector("list", nrow(starts))
  for (s in seq_len(nrow(starts))) {
    z <- as.numeric(starts[s, ])
    cur <- eval_at(z)
    mu <- 1e-3
    conv <- FALSE
    it <- 0L
    if (is.finite(cur$f)) {
      for (it in seq_len(maxit)) {
        Jo <- cur$jj$J[O, , drop = FALSE] * sw
        A <- crossprod(Jo) + diag(lambda, r)
        gneg <- as.numeric(crossprod(Jo, cur$res * sw)) - lambda * z
        improved <- FALSE
        while (mu < 1e12) {
          step <- tryCatch(solve(A + diag(mu * pmax(diag(A), 1e-12), r), gneg), error = function(e) NULL)
          if (!is.null(step)) {
            zn <- z + step
            nw <- eval_at(zn)
            if (is.finite(nw$f) && nw$f <= cur$f) {
              improved <- TRUE
              break
            }
          }
          mu <- mu * 4
        }
        if (!improved) {
          conv <- TRUE
          break
        }
        df <- cur$f - nw$f
        z <- zn
        cur <- nw
        mu <- max(mu / 4, 1e-15)
        if (sqrt(sum(step^2)) <= xtol * (1 + sqrt(sum(z^2))) || df <= 1e-15 * (1 + cur$f)) {
          conv <- TRUE
          break
        }
      }
    }
    sols[[s]] <- list(z = z, f = cur$f, iterations = it, converged = conv,
                      g = cur$jj$g, J = cur$jj$J)
  }
  fv <- vapply(sols, function(u) u$f, numeric(1))
  fv[!is.finite(fv)] <- Inf
  b <- which.min(fv)
  best <- sols[[b]]
  best$start <- b
  best$fvals <- fv
  best$spread <- NA_real_
  if (!is.null(tg)) {
    near <- which(fv <= fv[b] + 0.05 * (abs(fv[b]) + 1e-12))
    tv <- vapply(sols[near], function(u) u$g[tg[1]], numeric(1))
    best$spread <- max(tv) - min(tv)
  }
  best
}

# least-squares VAR(1) with intercept on the codes, eigenvalues outside the unit
# disk projected onto the unit circle
.ae_var1 <- function(S) {
  S <- as.matrix(S)
  K <- nrow(S)
  r <- ncol(S)
  Xd <- cbind(1, S[-K, , drop = FALSE])
  Yd <- S[-1, , drop = FALSE]
  B <- qr.solve(Xd, Yd)
  c0 <- as.numeric(B[1, ])
  A <- t(B[-1, , drop = FALSE])
  E <- Yd - Xd %*% B
  Sig <- crossprod(E) / max(1, K - 1 - (r + 1))
  ev <- eigen(A)
  rho <- max(Mod(ev$values))
  projected <- FALSE
  if (rho > 1) {
    lam <- ev$values / pmax(1, Mod(ev$values))
    A2 <- tryCatch(Re(ev$vectors %*% diag(lam, r) %*% solve(ev$vectors)), error = function(e) NULL)
    A <- if (is.null(A2)) A / rho else A2
    projected <- TRUE
  }
  list(c = c0, A = A, Sigma = Sig, rho = rho, projected = projected)
}

.ae_level_names <- function(level) paste0(format(100 * level, trim = TRUE), "%")

#' Forecasting with AESSA and VAESSA
#'
#' Point forecasts and (for [vaessa()] fits) predictive intervals.
#'
#' @section Methods:
#' \describe{
#'   \item{`"window"` (default)}{decoder inversion, one window per horizon.
#'   Horizon t <= L - r is predicted from the window that ends at N + t: its
#'   first L - t lags are observed, the code \eqn{\hat h} minimizes
#'   \eqn{\sum_{o} w_o (x_o - g_o(h))^2 + \lambda \|h\|^2} over the observed
#'   coordinates, and the forecast is the last lag of \eqn{g(\hat h)}. AESSA:
#'   \eqn{w_o = 1}, \eqn{\lambda = 0} (with a linear decoder, the minimum-norm
#'   recurrent forecast). VAESSA: \eqn{w_o = 1/\psi_o}, \eqn{\lambda = 1}, the
#'   MAP code; with a linear decoder this is the PrSSA window law of
#'   [prforecast()]. Horizons beyond L - r use the recursion below.}
#'   \item{`"recursive"`}{every horizon conditions on the last L - 1 values
#'   (observed values, then earlier forecasts) of its window.}
#'   \item{`"latent"`}{a VAR(1) with intercept is fitted to the codes
#'   \eqn{h_1, ..., h_K} (eigenvalues outside the unit disk are projected
#'   onto the unit circle), iterated from \eqn{h_K}, and decoded.}
#' }
#' The inversion is solved by Levenberg-Marquardt (analytic Jacobian of the
#' decoder) from several starts: the codes of the `starts` training windows
#' whose decoded observed coordinates fit best, the latent-VAR prediction, the
#' solution of the previous horizon and zero. The smallest objective wins;
#' `diagnostics$spread` is the range of the forecasts (first coordinate) over
#' starts whose objective is within 5 percent of the best (large values flag
#' phase-ambiguous inversions).
#'
#' @section Intervals (VAESSA only):
#' `"laplace"` (default for VAESSA): the code posterior is approximated by
#' \eqn{N(\hat h, (J_O' \Psi_O^{-1} J_O + I)^{-1})}, with \eqn{J_O} the
#' Jacobian of the observed coordinates at \eqn{\hat h}; the predictive
#' variance is \eqn{J_L \Sigma J_L' + \psi_L}. This is exact for a linear
#' decoder. For recursive horizons the earlier forecasts are treated as known,
#' which understates the uncertainty; use `"mc"` there. For the latent method
#' the VAR innovation covariance is propagated and mapped by the delta method.
#' `"mc"`: draws of the Laplace law are pushed through the decoder (window
#' horizons), recursive horizons are simulated path by path (each step
#' re-inverts the decoder; this is slow), and latent paths are simulated with
#' VAR innovations. Intervals are returned for scalar series (vector and matrix
#' input); for functional input only the point forecast is mapped to `funts`.
#'
#' @param fit an object from [aessa()] or [vaessa()].
#' @param h forecast horizon (number of steps).
#' @param method `"window"`, `"recursive"` or `"latent"` (see Section Methods).
#' @param lambda ridge weight of the inversion (`NULL`: 0 for AESSA, 1 for VAESSA).
#' @param interval `"none"`, `"laplace"` or `"mc"` (`NULL`: `"laplace"` for VAESSA, `"none"` for AESSA).
#' @param level confidence level(s) of the intervals.
#' @param npaths number of Monte Carlo draws or paths for `interval = "mc"`.
#' @param starts number of nearest training codes used as starting values.
#' @param maxit maximum number of Levenberg-Marquardt iterations per start.
#' @param seed seed for the Monte Carlo draws (`NULL`: current RNG state).
#' @return An object of class `"aeforecast"`:
#'   \item{mean}{point forecasts in the class of the input (vector or `ts`, h x p matrix, list, or `funts`)}
#'   \item{sd}{predictive standard deviations (vector or h x p matrix), or `NULL`}
#'   \item{lower, upper}{lists of interval bounds, one element per level (named as in [prforecast()]), or `NULL`}
#'   \item{latent}{h x r matrix of the codes behind each forecast}
#'   \item{method_used}{law used at each horizon}
#'   \item{diagnostics}{data frame: horizon, iterations, converged, objective, number of starts, spread}
#'   \item{var}{the latent VAR(1) (`c`, `A`, `Sigma`, spectral radius `rho`, `projected`)}
#'   \item{sd_working}{P x h predictive SDs in standardized working coordinates (all inputs)}
#'   \item{method, interval, level, lambda, h}{the settings}
#' @seealso [aessa()], [vaessa()], [prforecast()]
#' @examples
#' \donttest{
#' if (requireNamespace("torch", quietly = TRUE) && torch::torch_is_installed()) {
#'   set.seed(2)
#'   tt <- 1:200
#'   y <- tanh(3 * sin(2 * pi * tt / 25)) + rnorm(200, sd = 0.2)
#'   fit <- vaessa(y, L = 50, r = 2, restarts = 1, epochs = 200, seed = 1)
#'   fc <- aeforecast(fit, h = 12)
#'   fc
#'   aeforecast(fit, h = 12, method = "latent")$mean
#' }
#' }
#' @md
#' @export
aeforecast <- function(fit, h, method = c("window", "recursive", "latent"), lambda = NULL,
                       interval = NULL, level = 0.95, npaths = 500, starts = 3, maxit = 100,
                       seed = NULL) {
  .ae_need_torch()
  if (!inherits(fit, "aessa")) stop("`fit` must be an object from aessa() or vaessa().")
  method <- match.arg(method)
  vae <- inherits(fit, "vaessa")
  if (is.null(interval)) interval <- if (vae) "laplace" else "none"
  interval <- match.arg(interval, c("none", "laplace", "mc"))
  if (interval != "none" && !vae) {
    stop("Predictive intervals need a VAESSA fit (vaessa()); AESSA has no noise model.")
  }
  if (missing(h) || length(h) != 1 || !is.numeric(h) || h < 1 || h != round(h)) {
    stop("`h` must be a positive integer.")
  }
  if (!is.numeric(level) || any(level <= 0 | level >= 1)) stop("`level` must lie in (0, 1).")
  if (!is.numeric(npaths) || npaths < 10) stop("`npaths` must be at least 10.")
  h <- as.integer(h)
  npaths <- as.integer(npaths)
  if (is.null(lambda)) lambda <- if (vae) 1 else 0
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) stop("`lambda` must be a non-negative number.")
  if (!is.null(seed)) set.seed(seed)
  prep <- fit$prep
  spec <- fit$spec
  par <- fit$par
  L <- prep$L
  N <- prep$N
  K <- prep$K
  P <- prep$P
  r <- spec$r
  ri <- prep$ri
  rowmu <- prep$rowmu
  tg <- ri[, L]
  psi <- if (vae) fit$psi else rep(0, prep$D)
  wts <- if (vae) 1 / psi else rep(1, prep$D)
  S <- fit$scores
  G <- fit$records_hat
  var1 <- .ae_var1(S)
  Zext <- cbind(do.call(rbind, prep$Zs), matrix(NA_real_, P, h))
  mean_w <- matrix(NA_real_, P, h)
  sd_w <- matrix(NA_real_, P, h)
  Hs <- matrix(NA_real_, h, r)
  used <- character(h)
  dg <- data.frame(h = seq_len(h), iterations = NA_integer_, converged = NA, objective = NA_real_,
                   starts = NA_integer_, spread = NA_real_)
  draws <- if (interval == "mc") vector("list", h) else NULL
  Hw <- if (method == "window") min(h, L - r) else 0L
  if (method == "latent") {
    hc <- S[K, ]
    V <- matrix(0, r, r)
    for (t in seq_len(h)) {
      hc <- var1$c + as.numeric(var1$A %*% hc)
      V <- var1$A %*% V %*% t(var1$A) + var1$Sigma
      jj <- .ae_jac_R(par, hc, spec)
      mean_w[, t] <- jj$g[tg] + rowmu[tg]
      Hs[t, ] <- hc
      Jt <- jj$J[tg, , drop = FALSE]
      if (vae) sd_w[, t] <- sqrt(rowSums((Jt %*% V) * Jt) + psi[tg])
      used[t] <- "latent"
    }
    if (interval == "mc") {
      ch <- tryCatch(chol(var1$Sigma + diag(1e-12, r)), error = function(e) diag(sqrt(pmax(diag(var1$Sigma), 0)), r))
      Hp <- matrix(S[K, ], npaths, r, byrow = TRUE)
      for (t in seq_len(h)) {
        Hp <- sweep(Hp %*% t(var1$A), 2, var1$c, "+") + matrix(stats::rnorm(npaths * r), npaths, r) %*% ch
        Gp <- .ae_decode_R(par, Hp, spec)[, tg, drop = FALSE]
        draws[[t]] <- sweep(Gp, 2, rowmu[tg], "+") +
          matrix(stats::rnorm(npaths * P), npaths, P) * matrix(sqrt(psi[tg]), npaths, P, byrow = TRUE)
      }
    }
  } else {
    hv <- S[K, ]
    hprev <- NULL
    for (t in seq_len(h)) {
      Tt <- N + t
      nobs <- if (t <= Hw) L - t else L - 1L
      pp <- rep(seq_len(P), times = nobs)
      ll <- rep(seq_len(nobs), each = P)
      O <- ri[cbind(pp, ll)]
      xo <- Zext[cbind(pp, Tt - L + ll)] - rowmu[O]
      wo <- wts[O]
      hv <- var1$c + as.numeric(var1$A %*% hv)
      d2 <- as.numeric(((G[, O, drop = FALSE] - matrix(xo, K, length(O), byrow = TRUE))^2) %*% wo) +
        lambda * rowSums(S^2)
      st <- rbind(S[order(d2)[seq_len(min(starts, K))], , drop = FALSE], hv, hprev, rep(0, r))
      st <- unique(round(st, 12))
      inv <- .ae_invert(par, spec, xo, O, wo, lambda, st, tg = tg, maxit = maxit)
      mean_w[, t] <- inv$g[tg] + rowmu[tg]
      Zext[, Tt] <- mean_w[, t]
      Hs[t, ] <- inv$z
      hprev <- inv$z
      used[t] <- if (t <= Hw) "window" else "recursive"
      dg$iterations[t] <- inv$iterations
      dg$converged[t] <- inv$converged
      dg$objective[t] <- inv$f
      dg$starts[t] <- nrow(st)
      dg$spread[t] <- inv$spread
      if (vae) {
        Jo <- inv$J[O, , drop = FALSE]
        Sig <- tryCatch(solve(crossprod(Jo / sqrt(psi[O])) + diag(r)), error = function(e) diag(r))
        Jt <- inv$J[tg, , drop = FALSE]
        sd_w[, t] <- sqrt(rowSums((Jt %*% Sig) * Jt) + psi[tg])
        if (interval == "mc" && t <= Hw) {
          ch <- tryCatch(chol((Sig + t(Sig)) / 2), error = function(e) diag(sqrt(pmax(diag(Sig), 0)), r))
          Hp <- sweep(matrix(stats::rnorm(npaths * r), npaths, r) %*% ch, 2, inv$z, "+")
          Gp <- .ae_decode_R(par, Hp, spec)[, tg, drop = FALSE]
          draws[[t]] <- sweep(Gp, 2, rowmu[tg], "+") +
            matrix(stats::rnorm(npaths * P), npaths, P) * matrix(sqrt(psi[tg]), npaths, P, byrow = TRUE)
        }
      }
    }
    # Monte Carlo paths for the recursive horizons
    if (interval == "mc" && h > Hw) {
      Zobs <- do.call(rbind, prep$Zs)
      paths <- array(NA_real_, c(npaths, P, h))
      for (k in seq_len(npaths)) {
        Zk <- cbind(Zobs, matrix(NA_real_, P, h))
        if (Hw > 0) for (t in seq_len(Hw)) Zk[, N + t] <- draws[[t]][k, ]
        zp <- NULL
        for (t in seq.int(Hw + 1L, h)) {
          Tt <- N + t
          pp <- rep(seq_len(P), times = L - 1L)
          ll <- rep(seq_len(L - 1L), each = P)
          O <- ri[cbind(pp, ll)]
          xo <- Zk[cbind(pp, Tt - L + ll)] - rowmu[O]
          st <- rbind(Hs[t, ], zp)
          inv <- .ae_invert(par, spec, xo, O, wts[O], lambda, st, maxit = 50L)
          Jo <- inv$J[O, , drop = FALSE]
          Sig <- tryCatch(solve(crossprod(Jo / sqrt(psi[O])) + diag(r)), error = function(e) diag(r))
          ch <- tryCatch(chol((Sig + t(Sig)) / 2), error = function(e) diag(sqrt(pmax(diag(Sig), 0)), r))
          zz <- inv$z + as.numeric(stats::rnorm(r) %*% ch)
          val <- .ae_decode_R(par, matrix(zz, 1), spec)[1, tg] + rowmu[tg] + stats::rnorm(P) * sqrt(psi[tg])
          Zk[, Tt] <- val
          paths[k, , t] <- val
          zp <- inv$z
        }
      }
      for (t in seq.int(Hw + 1L, h)) draws[[t]] <- matrix(paths[, , t], npaths, P)
    }
  }
  # ---- map to the input class ---------------------------------------------------
  off <- c(0L, cumsum(prep$dims))
  to_blocks <- function(M) lapply(seq_along(prep$dims), function(b) M[off[b] + seq_len(prep$dims[b]), , drop = FALSE])
  idx <- N + seq_len(h)
  mean_out <- .ae_values_to_input(prep, to_blocks(mean_w), idx)
  sd_out <- lower <- upper <- NULL
  if (interval != "none") {
    sc_rows <- prep$scale[rep(seq_along(prep$dims), prep$dims)]
    sd_orig <- sd_w * sc_rows
    lw <- up <- vector("list", length(level))
    names(lw) <- names(up) <- .ae_level_names(level)
    for (i in seq_along(level)) {
      a <- (1 - level[i]) / 2
      if (interval == "laplace") {
        q <- stats::qnorm(1 - a)
        lo_w <- mean_w - q * sd_w
        up_w <- mean_w + q * sd_w
      } else {
        lo_w <- vapply(seq_len(h), function(t) apply(draws[[t]], 2, stats::quantile, probs = a, names = FALSE),
                       numeric(P))
        up_w <- vapply(seq_len(h), function(t) apply(draws[[t]], 2, stats::quantile, probs = 1 - a, names = FALSE),
                       numeric(P))
        lo_w <- matrix(lo_w, P, h)
        up_w <- matrix(up_w, P, h)
      }
      if (prep$scalar) {
        lw[[i]] <- .ae_values_to_input(prep, to_blocks(lo_w), idx)
        up[[i]] <- .ae_values_to_input(prep, to_blocks(up_w), idx)
      }
    }
    if (interval == "mc") {
      sd_w <- matrix(vapply(seq_len(h), function(t) apply(draws[[t]], 2, stats::sd), numeric(P)), P, h)
      sd_orig <- sd_w * sc_rows
    }
    if (prep$scalar) {
      sd_out <- if (prep$dat$type == "vector") as.numeric(sd_orig) else t(sd_orig)
      if (prep$dat$type == "matrix") colnames(sd_out) <- prep$dat$meta$colnames
      lower <- lw
      upper <- up
    }
  } else {
    sd_w <- NULL
  }
  out <- list(mean = mean_out, sd = sd_out, lower = lower, upper = upper, latent = Hs,
              method_used = used, diagnostics = if (method == "latent") NULL else dg,
              var = var1, sd_working = sd_w, method = method, interval = interval,
              level = level, lambda = lambda, h = h, model = if (vae) "VAESSA" else "AESSA")
  class(out) <- "aeforecast"
  out
}
