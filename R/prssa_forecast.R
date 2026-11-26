# =============================================================================
# PrSSA forecasting: window-conditional predictive law (h <= L - r), VAR(L-1)
# recursion with variance propagation (ridge LRF), and vector mode.
# =============================================================================

#' Forecasting for Probabilistic SSA
#'
#' Computes predictive distributions for the next \code{h} observations of the series
#' modelled by a \code{\link{prssa}} fit, with pointwise and simultaneous prediction bands.
#'
#' @details
#' \itemize{
#'   \item \code{"window"} (default): forecasting is imputation of the last lags of a window,
#'   with the Gaussian conditional law
#'   \eqn{N(\mu_F + W_F M_P^{-1}W_P^\top\Psi_P^{-1}(x_P - \mu_P),\ W_F M_P^{-1}W_F^\top + \Psi_F)},
#'   \eqn{M_P = I + W_P^\top\Psi_P^{-1}W_P}, of the future lags \eqn{F} given the observed
#'   lags \eqn{P}. By default (\code{joint = FALSE}) every horizon \eqn{t \le L - r} has its
#'   own window, the one ending at \eqn{N + t}: its last lag is predicted from its first
#'   \eqn{L - t} lags, the most observations a window can use (per-horizon marginal laws).
#'   With \code{joint = TRUE} all horizons \eqn{t \le H = \min(h, L - r)} come from the single
#'   window ending at \eqn{N + H}, conditioned on its first \eqn{L - H} lags; this is their
#'   joint law (returned in \code{joint_cov}, e.g. for simultaneous bands across horizons).
#'   The two coincide for \eqn{h = 1}, and horizon \eqn{t} of the marginal laws equals the
#'   last horizon of a joint call with \code{h = t}. Horizons beyond \eqn{L - r} use the
#'   recursive law.
#'   \item \code{"recursive"}: the one-step law (a ridge-regularized linear recurrence) is used
#'   as a VAR(\eqn{L-1}) working model; means are iterated and variances propagated through
#'   the impulse responses. As \eqn{\sigma^2 \to 0} the one-step mean tends to the SSA linear
#'   recurrence applied to the observed values (\code{Rssa::rforecast(base = "original")}).
#'   \item \code{"vector"}: SSA vector forecasting with the fitted loadings. The last
#'   (non-Hankelized) posterior-mean signal record \eqn{\hat s_K} is the seed, so the
#'   posterior shrinkage is applied once. It is then extended \eqn{h + L - 1} times by the
#'   noise-free V-forecast operator
#'   \eqn{P_V(y) = W (W_\nabla^\top\Psi_\nabla^{-1}W_\nabla)^{-1}W_\nabla^\top\Psi_\nabla^{-1}\bar y},
#'   where \eqn{\bar y} holds lags \eqn{2,\dots,L} of the current record and
#'   \eqn{W_\nabla}, \eqn{\Psi_\nabla} are the rows of lags \eqn{1,\dots,L-1}
#'   (the \eqn{\Psi^{-1}}-weighted least-squares projection onto \eqn{span(W)}), and the
#'   extended trajectory is diagonally averaged. With isotropic noise \eqn{P_V} is the
#'   V-forecast operator of the leading \eqn{r} eigenvectors, and as
#'   \eqn{\sigma^2 \to 0} the forecast is \code{Rssa::vforecast}. Only point forecasts are
#'   returned.
#' }
#' Pointwise bands are mean \eqn{\pm z_{(1+level)/2}} sd. Simultaneous bands use, for every
#' horizon, the \code{level} quantile of \eqn{\max_s |Z(s)|/sd(s)} over all reporting points
#' (curve grid points of all functional variables and all scalar series), estimated from
#' \code{nsim} Gaussian draws. For a single scalar series they coincide with pointwise bands.
#' Bands for functional variables are computed on the grid of the \code{funts} object and
#' also returned as \code{funts} objects (by least-squares projection on the basis).
#' With smoothing (\code{alpha > 0}) the raw-scale (second-stage) noise is used.
#'
#' @param fit An object of class \code{"prssa"}.
#' @param h Forecast horizon (positive integer).
#' @param method \code{"window"} (default), \code{"recursive"} or \code{"vector"}.
#' @param level Coverage levels of the prediction bands, in (0, 1).
#' @param bands Band types to compute: \code{"pointwise"} and/or \code{"simultaneous"} (both
#'   by default); \code{lower} and \code{upper} hold the first type requested.
#' @param nsim Number of Monte Carlo draws for simultaneous bands.
#' @param joint Window method only. \code{FALSE} (default): each horizon is conditioned on
#'   its own window (per-horizon marginal laws). \code{TRUE}: joint law of all window
#'   horizons from the largest window, with the cross-horizon covariance in \code{joint_cov}.
#'
#' @return An object of class \code{"prforecast"}, a list with
#' \item{mean}{Point forecasts, in the class of the input (numeric vector or \code{ts},
#'   \eqn{h \times p} matrix, list, or \code{\link{funts}} with \code{h} curves).}
#' \item{sd}{Predictive standard deviations at the reporting points: a vector, an
#'   \eqn{h \times p} matrix, or (functional variables) a grid \eqn{\times h} matrix per
#'   variable.}
#' \item{lower, upper}{Lists with one element per level (names such as \code{"95\%"}) in the
#'   class of the input, for the first band type in \code{bands}.}
#' \item{intervals}{All requested band types: \code{intervals$pointwise$lower} etc.}
#' \item{cov}{List over horizons of predictive covariance matrices in output coordinates
#'   (scalar values, basis coefficients, or grid values; blocks in input order).}
#' \item{quantile}{Multipliers used for the simultaneous bands (levels by horizons).}
#' \item{joint_cov}{With \code{method = "window"} and \code{joint = TRUE}: the joint predictive
#'   covariance of the window horizons \eqn{1,\dots,\min(h, L - r)} in output coordinates,
#'   horizon-major (all coordinates of horizon 1, then of horizon 2, ...); otherwise \code{NULL}.}
#' \item{method_used}{Law used at each horizon (\code{"window"}, \code{"recursive"} or
#'   \code{"vector"}).}
#' \item{level, bands, method, h, joint}{The arguments.}
#'
#' @seealso \code{\link{prssa}}, \code{\link{fforecast}}, \code{\link{fpredinterval}}
#'
#' @examples
#' set.seed(1)
#' tt <- 1:200
#' y <- 2 * sin(2 * pi * tt / 20) + cos(2 * pi * tt / 7) + rnorm(200, sd = 0.7)
#' fit <- prssa(y, L = 24, r = 4)
#' fc <- prforecast(fit, h = 10)
#' fc$mean
#' fc$lower[[2]]  # 95 percent pointwise band
#'
#' # functional time series: pointwise and simultaneous bands
#' u <- seq(0, 1, length.out = 25)
#' Xf <- outer(sin(pi * u), sin(2 * pi * tt[1:80] / 10)) + matrix(rnorm(25 * 80, sd = 0.2), 25)
#' Yf <- funts(Xf, fda::create.bspline.basis(c(0, 1), 9))
#' fcf <- prforecast(prssa(Yf, L = 10, r = 2), h = 2, nsim = 500)
#' fcf$quantile
#'
#' @export
prforecast <- function(fit, h, method = c("window", "recursive", "vector"), level = c(0.80, 0.95),
                       bands = c("pointwise", "simultaneous"), nsim = 1000, joint = FALSE) {
  if (!inherits(fit, "prssa")) stop("`fit` must be an object of class \"prssa\".")
  method <- match.arg(method)
  bands <- match.arg(bands, several.ok = TRUE)
  if (missing(h) || length(h) != 1 || h != round(h) || h < 1) stop("`h` must be a positive integer.")
  h <- as.integer(h)
  if (!is.numeric(level) || any(level <= 0 | level >= 1)) stop("`level` must lie in (0, 1).")
  if (!is.logical(joint) || length(joint) != 1 || is.na(joint)) stop("`joint` must be TRUE or FALSE.")
  dat <- fit$dat
  L <- fit$L
  dims <- dat$dims
  dtot <- sum(dims)
  N <- fit$N
  W <- fit$pred$W
  psi <- fit$pred$psi
  r <- ncol(W)
  ri <- .pr_record_rows(dims, L)
  Zc <- do.call(rbind, Map(function(b, m) b$Z - m, dat$blocks, fit$mean))
  mvec <- unlist(fit$mean, use.names = FALSE)
  means <- matrix(0, dtot, h)
  covs <- vector("list", h)
  used <- character(h)
  jcov <- NULL
  if (method == "vector") {
    means <- .pr_fc_vector(fit, W, psi, ri, h)
    used[] <- "vector"
  } else {
    hw <- if (method == "window") max(0L, min(h, L - r)) else 0L
    if (hw >= 1) {
      wc <- if (joint) .pr_fc_window(W, psi, ri, Zc, L, hw) else .pr_fc_window_marginal(W, psi, ri, Zc, L, hw)
      jcov <- wc$joint
      means[, seq_len(hw)] <- wc$mean
      covs[seq_len(hw)] <- wc$cov
      used[seq_len(hw)] <- "window"
    }
    if (h > hw) {
      rc <- .pr_fc_recursive(W, psi, ri, Zc, L, h)
      idx <- (hw + 1L):h
      means[, idx] <- rc$mean[, idx, drop = FALSE]
      covs[idx] <- rc$cov[idx]
      used[idx] <- "recursive"
    }
  }
  means <- means + mvec
  out <- .pr_fc_assemble(fit, means, covs, used, h, method, level, bands, nsim)
  if (!is.null(jcov)) {
    Tm <- .pr_fc_tmap(fit$dat)
    TT <- kronecker(diag(ncol(jcov) / nrow(Tm)), Tm)
    out$joint_cov <- TT %*% jcov %*% t(TT)
  }
  out$joint <- joint
  out
}

# window-conditional joint law of the next hw time vectors
.pr_fc_window <- function(W, psi, ri, Zc, L, hw) {
  N <- ncol(Zc)
  dtot <- nrow(ri)
  r <- ncol(W)
  np <- L - hw
  rowsP <- as.vector(ri[, seq_len(np), drop = FALSE])
  xP <- as.vector(Zc[, (N - np + 1L):N, drop = FALSE])
  WP <- W[rowsP, , drop = FALSE]
  WpP <- WP / psi[rowsP]
  MPinv <- solve(diag(r) + crossprod(WP, WpP))
  hP <- MPinv %*% crossprod(WpP, xP)
  mean <- matrix(0, dtot, hw)
  cov <- vector("list", hw)
  for (t in seq_len(hw)) {
    Ft <- ri[, np + t]
    WF <- W[Ft, , drop = FALSE]
    mean[, t] <- WF %*% hP
    cov[[t]] <- WF %*% MPinv %*% t(WF) + diag(psi[Ft], dtot)
  }
  Fall <- as.vector(ri[, np + seq_len(hw), drop = FALSE])
  WFa <- W[Fall, , drop = FALSE]
  joint <- WFa %*% MPinv %*% t(WFa) + diag(psi[Fall], length(Fall))
  list(mean = mean, cov = cov, joint = joint)
}

# per-horizon marginal window laws: horizon t conditions the window ending at
# N + t on its first L - t (observed) lags; equal to the last horizon of the
# joint law with hw = t
.pr_fc_window_marginal <- function(W, psi, ri, Zc, L, hw) {
  N <- ncol(Zc)
  dtot <- nrow(ri)
  r <- ncol(W)
  Ft <- ri[, L]
  WF <- W[Ft, , drop = FALSE]
  mean <- matrix(0, dtot, hw)
  cov <- vector("list", hw)
  for (t in seq_len(hw)) {
    np <- L - t
    rowsP <- as.vector(ri[, seq_len(np), drop = FALSE])
    xP <- as.vector(Zc[, (N - np + 1L):N, drop = FALSE])
    WP <- W[rowsP, , drop = FALSE]
    WpP <- WP / psi[rowsP]
    MPinv <- solve(diag(r) + crossprod(WP, WpP))
    mean[, t] <- WF %*% (MPinv %*% crossprod(WpP, xP))
    cov[[t]] <- WF %*% MPinv %*% t(WF) + diag(psi[Ft], dtot)
  }
  list(mean = mean, cov = cov, joint = NULL)
}

# one-step coefficients of the ridge LRF: y_{t+1} - m = W_L B (past L-1 values - m)
.pr_fc_onestep <- function(W, psi, ri, L) {
  r <- ncol(W)
  dtot <- nrow(ri)
  rowsP <- as.vector(ri[, seq_len(L - 1L), drop = FALSE])
  rowsL <- ri[, L]
  WP <- W[rowsP, , drop = FALSE]
  WpP <- WP / psi[rowsP]
  MPinv <- solve(diag(r) + crossprod(WP, WpP))
  WL <- W[rowsL, , drop = FALSE]
  list(B = MPinv %*% t(WpP), WL = WL, MPinv = MPinv, psiL = psi[rowsL],
       V1 = WL %*% MPinv %*% t(WL) + diag(psi[rowsL], dtot))
}

.pr_fc_recursive <- function(W, psi, ri, Zc, L, h) {
  N <- ncol(Zc)
  dtot <- nrow(ri)
  r <- ncol(W)
  os <- .pr_fc_onestep(W, psi, ri, L)
  path <- Zc[, (N - L + 2L):N, drop = FALSE]
  mean <- matrix(0, dtot, h)
  for (t in seq_len(h)) {
    np <- ncol(path)
    nw <- os$WL %*% (os$B %*% as.vector(path[, (np - L + 2L):np, drop = FALSE]))
    mean[, t] <- nw
    path <- cbind(path, nw)
  }
  Bl <- function(l) os$B[, (l - 1L) * dtot + seq_len(dtot), drop = FALSE]
  Cj <- lapply(seq_len(L - 1L), function(j) Bl(L - j) %*% os$WL)
  Phi <- vector("list", max(h - 1L, 0L))
  cov <- vector("list", h)
  cov[[1]] <- os$V1
  acc <- matrix(0, r, r)
  for (i in seq_len(h - 1L)) {
    Ph <- if (i <= L - 1L) Bl(L - i) else matrix(0, r, dtot)
    for (j in seq_len(min(i - 1L, L - 1L))) Ph <- Ph + Cj[[j]] %*% Phi[[i - j]]
    Phi[[i]] <- Ph
    PW <- Ph %*% os$WL
    acc <- acc + PW %*% os$MPinv %*% t(PW) + (Ph * rep(os$psiL, each = r)) %*% t(Ph)
    cov[[i + 1L]] <- os$V1 + os$WL %*% acc %*% t(os$WL)
  }
  list(mean = mean, cov = cov)
}

# vector mode: seed with the posterior-mean (shrunk once) signal records, then
# iterate the noise-free V-forecast operator of the fitted loadings,
#   P_V(y) = W (W_N' Psi_N^{-1} W_N)^{-1} W_N' Psi_N^{-1} y_D,
# i.e. the Psi^{-1}-weighted least-squares projection of the last L - 1 lags of
# the current record onto span(W) (rows N: lags 1..L-1, rows D: lags 2..L), and
# diagonally average the extended trajectory. No prior shrinkage enters the
# iteration, so the forecast is not damped step after step.
.pr_fc_vector <- function(fit, W, psi, ri, h) {
  L <- fit$L
  dims <- fit$dat$dims
  N <- fit$N
  sq <- sqrt(fit$psi)
  Shat <- (sq * fit$U) %*% ((fit$omega * fit$s) * t(fit$V))
  if (fit$alpha > 0) Shat <- .pr_apply_blockdiag(Shat, lapply(fit$dat$blocks, function(b) b$S), L, dims)
  rowsN <- as.vector(ri[, seq_len(L - 1L), drop = FALSE])
  rowsD <- as.vector(ri[, 2:L, drop = FALSE])
  WN <- W[rowsN, , drop = FALSE]
  WpN <- WN / psi[rowsN]
  G <- crossprod(WN, WpN)
  Ginv <- tryCatch(solve(G), error = function(e) .pr_pinv(G))
  PV <- W %*% (Ginv %*% t(WpN))
  ext <- matrix(0, nrow(Shat), h + L - 1L)
  last <- Shat[, ncol(Shat)]
  for (i in seq_len(h + L - 1L)) {
    last <- PV %*% last[rowsD]
    ext[, i] <- last
  }
  ser <- .pr_hankelize(cbind(Shat, ext), L, dims)
  do.call(rbind, lapply(ser, function(s) s[, N + seq_len(h), drop = FALSE]))
}

# convert grid/scalar values at the reporting points (list per block of
# n_b x h matrices) to the class of the input
.pr_eval_to_input <- function(dat, vals, idx) {
  conv_el <- function(el) {
    if (el$kind == "scalar") return(as.numeric(vals[[el$blocks]]))
    obj <- el$obj
    blks <- dat$blocks[el$blocks]
    vars <- vapply(blks, function(b) b$var, integer(1))
    Xl <- vector("list", length(blks))
    argv <- vector("list", length(blks))
    for (j in seq_along(blks)) {
      v <- vars[j]
      if (obj$dimSupp[[v]] == 1) {
        Xl[[j]] <- vals[[el$blocks[j]]]
        argv[[j]] <- as.numeric(obj$argval[[v]])
      } else {
        uu <- unique(obj$argval[[v]][, 1])
        vv <- unique(obj$argval[[v]][, 2])
        Xl[[j]] <- aperm(array(vals[[el$blocks[j]]], c(length(vv), length(uu), length(idx))), c(2, 1, 3))
        argv[[j]] <- list(uu, vv)
      }
    }
    tms <- .pr_funts_times(obj, idx)
    dn <- obj$dnames
    if (!is.null(dn) && is.list(dn)) dn <- dn[vars]
    vn <- obj$vnames
    if (!is.null(vn)) vn <- vn[vars]
    funts(X = Xl, basisobj = obj$basis[vars], argval = argv, method = "data",
          start = tms[1], end = tms[length(tms)], vnames = vn, dnames = dn, tname = obj$tname)
  }
  if (dat$type %in% c("vector", "matrix")) {
    return(.pr_as_input(dat, vals, idx))
  }
  out <- lapply(dat$meta$elements, conv_el)
  if (dat$type == "funts") return(out[[1]])
  names(out) <- vapply(dat$meta$elements, function(el) if (is.null(el$name)) "" else el$name, character(1))
  out
}

.pr_fc_assemble <- function(fit, means, covs, used, h, method, level, bands, nsim) {
  dat <- fit$dat
  idx <- fit$N + seq_len(h)
  ev <- .pr_eval_map(dat)
  E <- ev$E
  nb <- length(dat$blocks)
  split_blocks <- function(M, rows_list) lapply(rows_list, function(rr) M[rr, , drop = FALSE])
  wrows <- split(seq_len(sum(dat$dims)), rep(seq_len(nb), dat$dims))
  out <- list(mean = .pr_as_input(dat, split_blocks(means, wrows), idx))
  lev_names <- paste0(format(100 * level, trim = TRUE), "%")
  has_var <- method != "vector"
  meval <- E %*% means
  sdev <- matrix(NA_real_, nrow(E), h)
  qsim <- matrix(NA_real_, length(level), h, dimnames = list(lev_names, NULL))
  if (has_var) {
    for (t in seq_len(h)) {
      Vt <- covs[[t]]
      EV <- E %*% Vt
      sdev[, t] <- sqrt(pmax(rowSums(EV * E), 0))
      if ("simultaneous" %in% bands) {
        if (nrow(E) == 1) {
          qsim[, t] <- stats::qnorm(1 - (1 - level) / 2)
        } else {
          Rv <- E %*% .pr_msqrt(Vt)
          Zs <- Rv %*% matrix(stats::rnorm(ncol(Rv) * nsim), ncol(Rv))
          st <- apply(abs(Zs) / pmax(sdev[, t], .Machine$double.eps), 2, max)
          qsim[, t] <- stats::quantile(st, level, names = FALSE)
        }
      }
    }
  }
  sd_out <- .pr_fc_sd_output(dat, sdev, ev$rows)
  intervals <- list()
  for (bt in bands) {
    lo <- up <- list()
    for (k in seq_along(level)) {
      mult <- if (bt == "pointwise") {
        matrix(stats::qnorm(1 - (1 - level[k]) / 2), nrow(E), h)
      } else {
        matrix(qsim[k, ], nrow(E), h, byrow = TRUE)
      }
      lo[[lev_names[k]]] <- if (has_var) .pr_eval_to_input(dat, split_blocks(meval - mult * sdev, ev$rows), idx) else NULL
      up[[lev_names[k]]] <- if (has_var) .pr_eval_to_input(dat, split_blocks(meval + mult * sdev, ev$rows), idx) else NULL
    }
    intervals[[bt]] <- list(lower = lo, upper = up)
  }
  # covariances in output coordinates (scalar values, basis coefficients, grid values)
  Tmap <- .pr_fc_tmap(dat)
  out$sd <- sd_out
  out$lower <- intervals[[1]]$lower
  out$upper <- intervals[[1]]$upper
  out$intervals <- intervals
  out$cov <- if (has_var) lapply(covs, function(V) Tmap %*% V %*% t(Tmap)) else NULL
  out$quantile <- if (has_var && "simultaneous" %in% bands) qsim else NULL
  if (dat$type != "vector" && dat$type != "matrix") {
    out$grid <- list(mean = split_blocks(meval, ev$rows), sd = split_blocks(sdev, ev$rows))
  }
  out$method_used <- used
  out$level <- level
  out$bands <- bands
  out$method <- method
  out$h <- h
  class(out) <- "prforecast"
  out
}

# map from working to output coordinates (whitened coefficients -> basis coefficients)
.pr_fc_tmap <- function(dat) {
  Tmap <- matrix(0, sum(dat$dims), sum(dat$dims))
  off <- 0L
  for (b in seq_along(dat$blocks)) {
    blk <- dat$blocks[[b]]
    ii <- off + seq_len(blk$d)
    Tmap[ii, ii] <- if (blk$kind == "functional" && !isTRUE(blk$grid)) blk$Rinv else diag(blk$d)
    off <- off + blk$d
  }
  Tmap
}

.pr_fc_sd_output <- function(dat, sdev, rows) {
  if (dat$type == "vector") return(as.numeric(sdev[1, ]))
  if (dat$type == "matrix") {
    s <- t(sdev)
    colnames(s) <- dat$meta$colnames
    return(s)
  }
  per_block <- lapply(rows, function(rr) sdev[rr, , drop = FALSE])
  out <- lapply(dat$meta$elements, function(el) {
    if (el$kind == "scalar") return(as.numeric(per_block[[el$blocks]]))
    if (length(el$blocks) == 1) return(per_block[[el$blocks]])
    stats::setNames(per_block[el$blocks], dat$names[el$blocks])
  })
  if (dat$type == "funts") return(out[[1]])
  names(out) <- vapply(dat$meta$elements, function(el) if (is.null(el$name)) "" else el$name, character(1))
  out
}
