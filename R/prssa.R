# =============================================================================
# PrSSA / PrMSSA / PrFSSA / PrHSSA: fitting by closed form (isotropic noise)
# or ECME (block and heteroscedastic noise) of the window (composite)
# likelihood  sum_j log N(x_j; mu, W W' + Psi).
# =============================================================================

#' Probabilistic Singular Spectrum Analysis (PrSSA, PrMSSA, PrFSSA, PrHSSA)
#'
#' Fits the probabilistic SSA window model
#' \deqn{x_j = \mu + W h_j + \epsilon_j,\quad h_j \sim N(0, I_r),\quad \epsilon_j \sim N(0, \Psi),}
#' to the \eqn{L}-lagged records \eqn{x_j} (\eqn{j = 1, \dots, K = N - L + 1}) of a
#' univariate series (PrSSA), several series (PrMSSA, vertical stacking as in VMSSA/MFSSA),
#' a functional time series (PrFSSA) or a hybrid of these (PrHSSA), by maximizing the
#' composite (independence) log-likelihood over the overlapping windows.
#'
#' @details
#' \strong{Working space.} Scalar series enter as they are. Each functional variable of a
#' \code{\link{funts}} object is represented by its basis coefficients, whitened by the
#' symmetric square root of the exact \eqn{L^2} Gram matrix of its basis
#' (\code{fda::inprod}), so the Euclidean geometry of the working coordinates is the
#' \eqn{L^2} geometry of the curves. The noise level of PrFSSA is therefore a variance per
#' dimension of the span of the basis.
#'
#' \strong{Noise models.}
#' \itemize{
#'   \item \code{"iso"}: one variance \eqn{\sigma^2} for all coordinates. The maximizer has the
#'   closed form of Tipping and Bishop: the loadings span the leading \eqn{r} left singular
#'   vectors of the trajectory matrix and \eqn{\hat\sigma^2} is the mean discarded eigenvalue.
#'   With \code{mean = "zero"} this is exactly Basic SSA (FSSA, VMSSA for several series).
#'   \item \code{"block"}: one variance per series or functional variable (block). Fitted by
#'   ECME: the W-step is the truncated SVD of the noise-weighted trajectory matrix (weighted
#'   MSSA), the \eqn{\Psi}-step is an EM update. The log-likelihood never decreases.
#'   \item \code{"hetero"}: for functional variables, a noise variance that varies along the
#'   argument, \eqn{\log\sigma^2(s) = b(s)^\top\eta}, with a cubic B-spline basis \eqn{b} of
#'   dimension \code{hetero_df} and a second-difference (P-spline) penalty. The curves are
#'   represented by their values on the grid of the \code{funts} object. Scalar series in a
#'   hybrid list get one variance each.
#' }
#'
#' \strong{Mean.} \code{"plugin"} centres every series (or curve coordinate) at its sample
#' mean, so isotropic PrSSA is Basic SSA of the centred series; \code{"zero"} uses no
#' centring and reproduces uncentred Basic SSA exactly.
#'
#' \strong{Smoothing.} With \code{alpha > 0} (functional variables only), the model is fitted
#' to the half-smoothed records \eqn{S_\alpha x_j}, \eqn{S_\alpha = (I + \alpha Q)^{-1/2}} with
#' \eqn{Q} the second-derivative roughness operator (PrHSSA\eqn{_\alpha}). Its noise-free
#' limit is regularized HSSA. Loadings are reported as \eqn{S_\alpha \tilde W}. A second-stage
#' fit on the raw records, with loadings restricted to the smoothed signal subspace, supplies
#' the raw-scale noise used by \code{\link{prforecast}}.
#'
#' @param Y A numeric vector or \code{ts} (PrSSA); a numeric matrix with one series per
#'   column, or a list of numeric vectors (PrMSSA); a \code{\link{funts}} object, univariate
#'   or multivariate (PrFSSA); or a list mixing numeric vectors and \code{funts} objects
#'   (hybrid PrHSSA). All series must have the same length \eqn{N}.
#' @param L Window length, an integer with \eqn{2 \le L \le N - 1}.
#' @param r Number of latent factors (signal rank), \eqn{1 \le r < \min(D, K)} where \eqn{D}
#'   is the dimension of a record.
#' @param noise Noise model: \code{"iso"} (default), \code{"block"} or \code{"hetero"}.
#' @param mean Centring: \code{"plugin"} (default, sample mean of each series) or
#'   \code{"zero"} (no centring; exact Basic SSA).
#' @param alpha Non-negative smoothing parameter of the half-smoothing operator
#'   (functional variables only). \code{0} (default) means no smoothing.
#' @param hetero_df Number of basis functions for \eqn{\log\sigma^2(s)} when
#'   \code{noise = "hetero"}: 1 gives a constant (isotropic on the grid), 2 or 3 a polynomial,
#'   4 or more a cubic B-spline basis.
#' @param hetero_lambda Weight of the second-difference penalty on \eqn{\eta}, in units of the
#'   composite-likelihood weight \eqn{KL/2} of one grid point. \code{NULL} (default) uses 1;
#'   \code{0} gives an unpenalized regression spline.
#' @param control A list of control parameters: \code{maxit} (default 500) and \code{tol}
#'   (default 1e-8, relative change of the noise variances) for ECME; \code{trace}
#'   (\code{FALSE}; if \code{TRUE} the composite log-likelihood after every half-step is kept
#'   in \code{$trace}); \code{floor} (1e-6, lower bound of a variance relative to the average
#'   variance of its block); \code{nextra} (10, extra eigenvalues kept for the scree plot);
#'   \code{svd} (\code{"auto"}, \code{"base"} or \code{"svds"}).
#'
#' @return An object of class \code{"prssa"}, a list with components
#' \item{W}{\eqn{D \times r} loading matrix in working coordinates (records stacked
#'   series/variable-outer, coordinate next, lag inner).}
#' \item{sigma2}{Estimated noise variance(s): a number (\code{"iso"}), a named vector with
#'   one value per block (\code{"block"}), or a named list whose functional elements are
#'   \eqn{\sigma^2(s)} on the grid (\code{"hetero"}).}
#' \item{psi}{The diagonal of \eqn{\Psi} (length \eqn{D}).}
#' \item{mu}{The mean record (length \eqn{D}); \code{mean} holds the mean of each series in
#'   working coordinates.}
#' \item{eigenvalues}{Leading eigenvalues of the sample record covariance (\code{"iso"}), or
#'   of its noise-whitened version \eqn{\Psi^{-1/2}\hat\Sigma\Psi^{-1/2}} (other noise
#'   models).}
#' \item{lambda, omega}{Noise-whitened eigenvalues \eqn{\tilde\lambda_k} of the retained
#'   components and posterior shrinkage factors \eqn{\omega_k = 1 - 1/\tilde\lambda_k}.}
#' \item{loglik}{Maximized composite log-likelihood (penalized for \code{"hetero"} with a
#'   penalty).}
#' \item{df}{Nominal number of free parameters \eqn{Dr - r(r-1)/2 + } (number of noise
#'   parameters).}
#' \item{iterations, converged}{ECME iterations and convergence flag (0 and \code{TRUE} for
#'   the closed form).}
#' \item{call}{The matched call.}
#' and further components used by \code{\link{prreconstruct}}, \code{\link{prforecast}} and
#' \code{\link{prrank}} (\code{r}, \code{L}, \code{K}, \code{N}, \code{D}, the singular
#' triples \code{U}, \code{s}, \code{V}, the prediction parameters \code{pred} and the
#' prepared data).
#'
#' @references
#' \enumerate{
#'   \item Tipping, M. E., & Bishop, C. M. (1999). Probabilistic principal component analysis.
#'   \emph{Journal of the Royal Statistical Society: Series B}, \strong{61}(3), 611-622.
#'   \item Haghbin, H., Najibi, S. M., Mahmoudvand, R., Trinka, J., & Maadooliat, M. (2021).
#'   Functional singular spectrum analysis. \emph{Stat}, \strong{10}(1), e330.
#' }
#' @seealso \code{\link{prreconstruct}}, \code{\link{prforecast}}, \code{\link{prrank}},
#'   \code{\link{fssa}}
#'
#' @examples
#' set.seed(1)
#' N <- 200
#' tt <- 1:N
#' y <- 2 * sin(2 * pi * tt / 20) + cos(2 * pi * tt / 7) + rnorm(N, sd = 0.7)
#' fit <- prssa(y, L = 24, r = 4)
#' fit
#'
#' # PrMSSA with one noise level per series
#' y2 <- 1.5 * sin(2 * pi * tt / 20 + 0.5) + rnorm(N, sd = 2)
#' fit2 <- prssa(cbind(y1 = y, y2 = y2), L = 24, r = 2, noise = "block")
#' fit2$sigma2
#'
#' # PrFSSA of a functional time series on a B-spline basis
#' u <- seq(0, 1, length.out = 25)
#' Xf <- outer(sin(pi * u), sin(2 * pi * tt[1:80] / 10)) + matrix(rnorm(25 * 80, sd = 0.2), 25)
#' Yf <- funts(Xf, fda::create.bspline.basis(c(0, 1), 9))
#' fitf <- prssa(Yf, L = 10, r = 2)
#' fitf$sigma2
#'
#' @export
prssa <- function(Y, L, r, noise = c("iso", "block", "hetero"), mean = c("plugin", "zero"),
                  alpha = 0, hetero_df = 8, hetero_lambda = NULL,
                  control = list(maxit = 500, tol = 1e-8)) {
  cl <- match.call()
  noise <- match.arg(noise)
  mean <- match.arg(mean)
  ctrl <- .pr_control(control)
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha < 0) stop("`alpha` must be a non-negative number.")
  dat <- .pr_prepare(Y, noise, alpha)
  N <- dat$N
  if (missing(L) || length(L) != 1 || L != round(L) || L < 2 || L > N - 1) {
    stop("`L` must be an integer with 2 <= L <= N - 1.")
  }
  L <- as.integer(L)
  K <- N - L + 1L
  dims <- dat$dims
  D <- sum(dims) * L
  if (missing(r) || length(r) != 1 || r != round(r) || r < 1 || r >= min(D, K)) {
    stop("`r` must be an integer with 1 <= r < min(D, K).")
  }
  r <- as.integer(r)
  if (alpha > 0 && !any(vapply(dat$blocks, function(b) !is.null(b$S), logical(1)))) {
    stop("Smoothing (alpha > 0) needs at least one functional (`funts`) variable.")
  }
  # mean of each series in working coordinates and the centred trajectory
  mlist <- lapply(dat$blocks, function(b) if (mean == "plugin") rowMeans(b$Z) else rep(0, b$d))
  Zc <- Map(function(b, m) b$Z - m, dat$blocks, mlist)
  X <- .pr_trajectory(Zc, L)
  Slist <- lapply(dat$blocks, function(b) b$S)
  smooth <- alpha > 0
  Xs <- if (smooth) .pr_apply_blockdiag(X, Slist, L, dims) else X
  nz <- .pr_noise_structure(dat, noise, L, K, hetero_df, hetero_lambda)
  if (noise != "iso" && any(vapply(nz$groups, function(g) {
    length(which(nz$unit %in% g$units))
  }, numeric(1)) <= r)) {
    warning("Some noise blocks have dimension <= r; their noise levels may not be identified.")
  }
  est <- if (nz$iso) .pr_fit_iso(Xs, r, ctrl) else .pr_fit_ecme(Xs, r, nz, ctrl)
  if (ncol(est$W) == 0) stop("No component exceeds the noise level; reduce `r` or check the data.")
  rr <- ncol(est$W)
  if (rr < r) warning(sprintf("Only %d of the %d requested components exceed the noise level.", rr, r))
  ll <- .pr_cl_loglik(Xs, est$W, est$psi) - .pr_penalty(nz, est$par)
  # loadings in the reported metric and the prediction (raw-scale) parameters
  W <- if (smooth) .pr_apply_blockdiag(est$W, Slist, L, dims) else est$W
  pred <- list(W = est$W, psi = est$psi)
  stage2 <- NULL
  if (smooth) {
    Fb <- .pr_orth(.pr_apply_blockdiag(sqrt(est$psi) * est$U, Slist, L, dims))
    s2 <- if (nz$iso) .pr_fit_iso_restricted(X, Fb) else .pr_fit_ecme(X, rr, nz, ctrl, Fb = Fb)
    pred <- list(W = s2$W, psi = s2$psi)
    stage2 <- list(W = s2$W, psi = s2$psi, sigma2 = .pr_sigma2_report(s2$unitvar, nz, dat),
                   loglik = .pr_cl_loglik(X, s2$W, s2$psi) - .pr_penalty(nz, s2$par))
  }
  ntot <- .pr_noise_npar(nz)
  fit <- list(
    W = W, sigma2 = .pr_sigma2_report(est$unitvar, nz, dat), psi = est$psi,
    mu = unlist(lapply(seq_along(mlist), function(b) rep(mlist[[b]], each = L))),
    mean = stats::setNames(mlist, dat$names),
    eigenvalues = est$eig, eig_scale = est$eig_scale,
    lambda = est$lambda, omega = 1 - 1 / est$lambda,
    loglik = ll, df = D * rr - rr * (rr - 1) / 2 + ntot, npar_noise = ntot,
    iterations = est$iterations, converged = est$converged, trace = est$trace,
    U = est$U, s = est$s, V = est$V, pred = pred, stage2 = stage2,
    noise_par = est$par, r = rr, r_requested = r, L = L, K = K, N = N, D = D,
    noise = noise, mean_type = mean, alpha = alpha, hetero_df = hetero_df,
    dat = dat, nz = nz, ctrl = ctrl, call = cl
  )
  class(fit) <- "prssa"
  fit
}

# ------------------------------------------------------------------------------
# Composite log-likelihood (Woodbury form; never forms a D x D matrix)
# ------------------------------------------------------------------------------
.pr_cl_loglik <- function(X, W, psi) {
  D <- nrow(X)
  K <- ncol(X)
  r <- ncol(W)
  q0 <- sum(X^2 / psi)
  if (r == 0) return(-0.5 * (K * D * log(2 * pi) + K * sum(log(psi)) + q0))
  Wp <- W / psi
  M <- diag(r) + crossprod(W, Wp)
  A <- crossprod(Wp, X)
  ld <- sum(log(psi)) + as.numeric(determinant(M, logarithm = TRUE)$modulus)
  -0.5 * (K * D * log(2 * pi) + K * ld + q0 - sum(A * solve(M, A)))
}

.pr_penalty <- function(nz, par) {
  pen <- 0
  for (g in seq_along(nz$groups)) {
    G <- nz$groups[[g]]
    if (G$type == "loglin" && G$lambda > 0) {
      eta <- par[[g]]
      pen <- pen + 0.5 * G$lambda * sum(eta * (G$P %*% eta))
    }
  }
  pen
}

.pr_sigma2_report <- function(unitvar, nz, dat) {
  if (nz$iso) return(unitvar[1])
  if (all(vapply(nz$groups, function(g) g$type == "const", logical(1)))) {
    return(stats::setNames(vapply(nz$groups, function(g) unitvar[g$units], numeric(1)),
                           vapply(nz$groups, function(g) g$name, character(1))))
  }
  out <- lapply(nz$groups, function(g) {
    if (g$type == "const") return(unitvar[g$units])
    v <- unitvar[g$units]
    names(v) <- format(dat$blocks[[g$block]]$s, digits = 4)
    v
  })
  names(out) <- vapply(nz$groups, function(g) g$name, character(1))
  out
}

# ------------------------------------------------------------------------------
# Isotropic closed form (Tipping-Bishop on the composite likelihood)
# ------------------------------------------------------------------------------
.pr_fit_iso <- function(X, r, ctrl) {
  D <- nrow(X)
  K <- ncol(X)
  k <- min(r + ctrl$nextra, min(D, K))
  sv <- .pr_svd(X, k, ctrl)
  kap <- sv$d^2 / K
  sig2 <- (sum(X^2) / K - sum(kap[seq_len(r)])) / (D - r)
  if (!(sig2 > 0)) stop("The estimated noise variance is not positive; reduce `r`.")
  keep <- which(kap[seq_len(r)] > sig2)
  U <- sv$u[, keep, drop = FALSE]
  W <- sweep(U, 2, sqrt(kap[keep] - sig2), "*")
  list(W = W, psi = rep(sig2, D), U = U, V = sv$v[, keep, drop = FALSE],
       s = sv$d[keep] / sqrt(sig2), lambda = kap[keep] / sig2, eig = kap, eig_scale = sig2,
       iterations = 0L, converged = TRUE, trace = NULL, unitvar = sig2, par = list(log(sig2)))
}

# isotropic raw-scale fit with loadings restricted to span(Fb) (stage 2 of smoothing)
.pr_fit_iso_restricted <- function(X, Fb) {
  K <- ncol(X)
  D <- nrow(X)
  rf <- ncol(Fb)
  A <- crossprod(Fb, X)
  sig2 <- (sum(X^2) - sum(A^2)) / (K * (D - rf))
  e <- eigen(tcrossprod(A) / K, symmetric = TRUE)
  keep <- which(e$values > sig2)
  U <- Fb %*% e$vectors[, keep, drop = FALSE]
  list(W = sweep(U, 2, sqrt(e$values[keep] - sig2), "*"), psi = rep(sig2, D),
       unitvar = sig2, par = list(log(sig2)))
}

# ------------------------------------------------------------------------------
# ECME for block / heteroscedastic noise
# ------------------------------------------------------------------------------
# W-step: exact maximizer over W for fixed diagonal Psi (optionally with span(W) in span(Fb))
.pr_wstep <- function(X, r, psi, ctrl, Fb = NULL, k_extra = 0) {
  K <- ncol(X)
  sq <- sqrt(psi)
  Xt <- X / sq
  if (is.null(Fb)) {
    sv <- .pr_svd(Xt, min(r + k_extra, min(dim(Xt))), ctrl)
    lam_all <- sv$d^2 / K
    lam <- lam_all[seq_len(r)]
    keep <- which(lam > 1)
    U <- sv$u[, keep, drop = FALSE]
    V <- sv$v[, keep, drop = FALSE]
    s <- sv$d[keep]
  } else {
    Q <- .pr_orth(Fb / sq)
    A <- crossprod(Q, Xt)
    e <- eigen(tcrossprod(A) / K, symmetric = TRUE)
    lam_all <- e$values
    lam <- lam_all[seq_len(min(r, length(lam_all)))]
    keep <- which(lam > 1)
    U <- Q %*% e$vectors[, keep, drop = FALSE]
    s <- sqrt(K * lam[keep])
    V <- sweep(crossprod(Xt, U), 2, s, "/")
  }
  W <- sq * sweep(U, 2, sqrt(lam[keep] - 1), "*")
  list(W = W, U = U, V = V, s = s, lambda = lam[keep], eig = lam_all)
}

# diagonal of the EM residual matrix E[(x - W h)(x - W h)'] averaged over windows
.pr_estep_diag <- function(X, W, psi, rs2) {
  K <- ncol(X)
  r <- ncol(W)
  if (r == 0) return(rs2)
  Wp <- W / psi
  Minv <- solve(diag(r) + crossprod(W, Wp))
  BX <- Minv %*% crossprod(Wp, X)
  XB <- X %*% t(BX) / K
  Ehh <- Minv + tcrossprod(BX) / K
  rs2 - 2 * rowSums(W * XB) + rowSums((W %*% Ehh) * W)
}

.pr_init_noise <- function(X, r, nz, ctrl, rows_by_unit) {
  K <- ncol(X)
  rs <- rowSums(X^2)
  unitvar <- numeric(length(nz$unit_n))
  par <- vector("list", length(nz$groups))
  floors <- numeric(length(nz$groups))
  blockwise <- function(rows) {
    Xg <- X[rows, , drop = FALSE]
    tot <- sum(rs[rows])
    Dg <- length(rows)
    v <- if (Dg > r) {
      d <- svd(Xg, nu = 0, nv = 0)$d[seq_len(r)]
      (tot - sum(d^2)) / (K * (Dg - r))
    } else {
      tot / (K * Dg) / 2
    }
    c(v, tot / (K * Dg))
  }
  for (g in seq_along(nz$groups)) {
    G <- nz$groups[[g]]
    if (G$type == "const") {
      rows <- unlist(rows_by_unit[G$units], use.names = FALSE)
      vb <- blockwise(rows)
      floors[g] <- ctrl$floor * vb[2]
      v <- max(vb[1], floors[g])
      par[[g]] <- log(v)
      unitvar[G$units] <- v
    } else {
      vb <- vapply(G$units, function(u) blockwise(rows_by_unit[[u]]), numeric(2))
      floors[g] <- ctrl$floor * mean(vb[2, ])
      s20 <- pmax(vb[1, ], floors[g])
      eta <- qr.solve(G$Bd, log(s20))
      par[[g]] <- eta
      unitvar[G$units] <- exp(drop(G$Bd %*% eta))
    }
  }
  list(par = par, unitvar = unitvar, floors = floors)
}

# penalized Newton step for log sigma^2 = B eta (concave gamma-type regression)
.pr_loglin_update <- function(eta, ebar, n, Bd, P, lambda, K) {
  Qf <- function(et) {
    lp <- drop(Bd %*% et)
    -(K / 2) * sum(n * (lp + ebar * exp(-lp))) - 0.5 * lambda * sum(et * (P %*% et))
  }
  for (it in seq_len(100)) {
    lp <- drop(Bd %*% eta)
    u <- ebar * exp(-lp)
    gr <- -(K / 2) * colSums(Bd * (n * (1 - u))) - lambda * drop(P %*% eta)
    H <- -(K / 2) * crossprod(Bd, Bd * (n * u)) - lambda * P
    st <- -solve(H, gr)
    tt <- 1
    q0 <- Qf(eta)
    while (Qf(eta + tt * st) < q0 && tt > 1e-10) tt <- tt / 2
    eta <- eta + tt * st
    if (max(abs(tt * st)) < 1e-12) break
  }
  eta
}

.pr_psistep <- function(e, nz, par, K, floors) {
  esum <- as.numeric(rowsum(e, nz$unit, reorder = TRUE))
  ebar <- esum / nz$unit_n
  unitvar <- numeric(length(ebar))
  for (g in seq_along(nz$groups)) {
    G <- nz$groups[[g]]
    if (G$type == "const") {
      v <- max(sum(esum[G$units]) / sum(nz$unit_n[G$units]), floors[g])
      par[[g]] <- log(v)
      unitvar[G$units] <- v
    } else {
      eb <- pmax(ebar[G$units], floors[g])
      eta <- .pr_loglin_update(par[[g]], eb, nz$unit_n[G$units], G$Bd, G$P, G$lambda, K)
      par[[g]] <- eta
      unitvar[G$units] <- exp(drop(G$Bd %*% eta))
    }
  }
  list(par = par, unitvar = unitvar)
}

.pr_fit_ecme <- function(X, r, nz, ctrl, Fb = NULL) {
  K <- ncol(X)
  D <- nrow(X)
  rows_by_unit <- split(seq_len(D), nz$unit)
  init <- .pr_init_noise(X, r, nz, ctrl, rows_by_unit)
  par <- init$par
  unitvar <- init$unitvar
  rs2 <- rowSums(X^2) / K
  tr <- if (isTRUE(ctrl$trace)) numeric(0) else NULL
  conv <- FALSE
  it <- 0L
  for (it in seq_len(ctrl$maxit)) {
    psi <- unitvar[nz$unit]
    ws <- .pr_wstep(X, r, psi, ctrl, Fb)
    if (!is.null(tr)) tr <- c(tr, .pr_cl_loglik(X, ws$W, psi) - .pr_penalty(nz, par))
    e <- .pr_estep_diag(X, ws$W, psi, rs2)
    new <- .pr_psistep(e, nz, par, K, init$floors)
    if (!is.null(tr)) {
      tr <- c(tr, .pr_cl_loglik(X, ws$W, new$unitvar[nz$unit]) - .pr_penalty(nz, new$par))
    }
    delta <- max(abs(new$unitvar - unitvar) / unitvar)
    par <- new$par
    unitvar <- new$unitvar
    if (delta < ctrl$tol) {
      conv <- TRUE
      break
    }
  }
  psi <- unitvar[nz$unit]
  ws <- .pr_wstep(X, r, psi, ctrl, Fb, k_extra = if (is.null(Fb)) ctrl$nextra else 0)
  list(W = ws$W, psi = psi, U = ws$U, V = ws$V, s = ws$s, lambda = ws$lambda, eig = ws$eig,
       eig_scale = 1, iterations = it, converged = conv, trace = tr, unitvar = unitvar, par = par)
}

# rebuild the (centred, smoothed when alpha > 0) working trajectory of a fit
.pr_fit_trajectory <- function(fit, smoothed = TRUE) {
  Zc <- Map(function(b, m) b$Z - m, fit$dat$blocks, fit$mean)
  X <- .pr_trajectory(Zc, fit$L)
  if (smoothed && fit$alpha > 0) {
    X <- .pr_apply_blockdiag(X, lapply(fit$dat$blocks, function(b) b$S), fit$L, fit$dat$dims)
  }
  X
}

# the sigma -> 0 limit of a fit with the signal subspace held fixed: noise
# scaled by tau, loadings  Psi^{1/2} U (Lambda - tau I)^{1/2}  (used in tests
# of the reductions to SSA forecasting)
.pr_noise_limit <- function(fit, tau) {
  sq <- sqrt(fit$psi)
  lim <- fit
  lim$pred$psi <- tau * fit$psi
  lim$pred$W <- sq * sweep(fit$U, 2, sqrt(fit$lambda - tau), "*")
  lim$psi <- lim$pred$psi
  lim$lambda <- fit$lambda / tau
  lim$omega <- 1 - 1 / lim$lambda
  lim$s <- fit$s / sqrt(tau)
  lim$U <- fit$U
  lim
}
