# ------------------------------------------------------------------------------
# Internal helpers of the AESSA / VAESSA family (autoencoder SSA).
#
# Conventions
#   * Records (lag vectors) are stored as ROWS of a K x D matrix; the record of
#     coordinate p (blocks stacked, coordinates of a block, lag inner) at lag l
#     is column .pr_record_rows(dims, L)[p, l] (the PrSSA layout).
#   * Working values are centred by the mean convention, divided by a scale per
#     block and (convention M2 only) centred again by the mean record.
#   * The network parameters are kept as a flat named list of base-R arrays
#     (`par`); torch tensors exist only during training. Decoding, Jacobians and
#     decoder inversion run in base R (double precision) on `par`.
# ------------------------------------------------------------------------------

.ae_need_torch <- function() {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("AESSA/VAESSA need the 'torch' package. Install it with install.packages(\"torch\") ",
         "and then run torch::install_torch().", call. = FALSE)
  }
  if (!isTRUE(torch::torch_is_installed())) {
    stop("The 'torch' package is installed but LibTorch is not: run torch::install_torch().",
         call. = FALSE)
  }
  invisible(TRUE)
}

# activation functions: base-R value, derivative (from pre-activation and value), torch value
.ae_act <- function(name) {
  switch(name,
    tanh = list(R = tanh, dR = function(pre, a) 1 - a^2,
                T = function(x) torch::torch_tanh(x)),
    sigmoid = list(R = function(x) 1 / (1 + exp(-x)), dR = function(pre, a) a * (1 - a),
                   T = function(x) torch::torch_sigmoid(x)),
    softplus = list(R = function(x) ifelse(x > 30, x, log1p(exp(pmin(x, 30)))),
                    dR = function(pre, a) 1 / (1 + exp(-pre)),
                    T = function(x) torch::nnf_softplus(x)),
    stop("Unknown activation: ", name)
  )
}

# ------------------------------------------------------------------------------
# Data preparation (re-uses the PrSSA input layer: vector/ts, matrix, funts, list)
# ------------------------------------------------------------------------------
.ae_prepare <- function(Y, L, mean_type, scale) {
  dat <- .pr_prepare(Y, "iso", 0)
  N <- dat$N
  if (missing(L) || length(L) != 1 || !is.numeric(L) || L != round(L) || L < 2 || L > N - 1) {
    stop("`L` must be an integer with 2 <= L <= N - 1.")
  }
  L <- as.integer(L)
  K <- N - L + 1L
  dims <- dat$dims
  D <- sum(dims) * L
  mlist <- lapply(dat$blocks, function(b) if (mean_type == "M1") rowMeans(b$Z) else rep(0, b$d))
  Zc <- Map(function(b, m) b$Z - m, dat$blocks, mlist)
  sc <- if (isTRUE(scale)) vapply(Zc, function(z) sqrt(mean(z^2)), numeric(1)) else rep(1, length(Zc))
  sc[!is.finite(sc) | sc <= 0] <- 1
  Zs <- Map(function(z, s) z / s, Zc, sc)
  Xr <- t(.pr_trajectory(Zs, L))                       # K x D, records as rows
  rowmu <- if (mean_type == "M2") colMeans(Xr) else rep(0, D)
  Xr <- sweep(Xr, 2, rowmu)
  ri <- .pr_record_rows(dims, L)
  blk <- rep(seq_along(dims), dims * L)
  list(dat = dat, N = N, L = L, K = K, D = D, dims = dims, P = sum(dims),
       mlist = mlist, scale = sc, Zs = Zs, X = Xr, rowmu = rowmu, ri = ri, blk = blk,
       mean_type = mean_type, scalar = all(vapply(dat$blocks, function(b) b$kind == "scalar", logical(1))))
}

# noise units of VAESSA: one (iso) or one per block (block)
.ae_units <- function(prep, noise) {
  unit <- if (noise == "iso") rep(1L, prep$D) else prep$blk
  nu <- max(unit)
  A <- matrix(0, prep$D, nu)
  A[cbind(seq_len(prep$D), unit)] <- 1
  nm <- if (noise == "iso") "all" else prep$dat$names
  list(unit = unit, n = nu, A = A, names = nm)
}

# standardised working records (rows) -> series of the input class (adds scale and mean)
# Tm: n x D matrix of centred standardised records (without rowmu); idx: time indices
.ae_records_to_input <- function(prep, Tm, idx, add_mean = TRUE) {
  Tm <- if (add_mean) sweep(Tm, 2, prep$rowmu, "+") else Tm
  Tm <- sweep(Tm, 2, prep$scale[prep$blk], "*")
  ser <- .pr_hankelize(t(Tm), prep$L, prep$dims)
  if (add_mean) ser <- Map(function(s, m) s + m, ser, prep$mlist)
  .pr_as_input(prep$dat, ser, idx)
}

# the mean part (constant series of the plug-in mean, plus the Hankelised mean
# record under convention M2) in the input class
.ae_mean_input <- function(prep) {
  K <- prep$K
  Tm <- matrix(0, K, prep$D)
  .ae_records_to_input(prep, Tm, seq_len(prep$N), add_mean = TRUE)
}

# series values (P x n working coordinates per block, standardised) -> original
# scale and class; vals: list over blocks of d_b x n matrices in standardised
# centred units, rows of rowmu already added by the caller
.ae_values_to_input <- function(prep, vals, idx) {
  ser <- Map(function(v, s, m) v * s + m, vals, as.list(prep$scale), prep$mlist)
  .pr_as_input(prep$dat, ser, idx)
}

# ------------------------------------------------------------------------------
# Linear solutions used as starting values
# ------------------------------------------------------------------------------
# Basic SSA of the records (rows of X): W = U_r S_r / sqrt(n), E = sqrt(n) S_r^{-1} U_r'
.ae_linear_ssa <- function(X, r) {
  n <- nrow(X)
  sv <- svd(X, nu = 0, nv = r)
  d <- sv$d[seq_len(r)]
  if (any(d <= 1e-12 * max(sv$d))) stop("The trajectory matrix has rank < r; reduce `r`.")
  U <- sv$v[, seq_len(r), drop = FALSE]
  list(W = sweep(U, 2, d / sqrt(n), "*"), E = t(sweep(U, 2, sqrt(n) / d, "*")),
       kappa = sv$d^2 / n)
}

# probabilistic PCA (PrSSA, iso or block noise) of the records; returns loadings
# rotated so that W' Psi^-1 W is diagonal, the exact (diagonal) posterior
# encoder and the unit noise variances
.ae_linear_ppca <- function(X, r, units, dat, L, noise) {
  n <- nrow(X)
  D <- ncol(X)
  if (noise == "iso") {
    sv <- svd(X, nu = 0, nv = min(r, D))
    kap <- sv$d^2 / n
    s2 <- (sum(X^2) / n - sum(kap[seq_len(r)])) / (D - r)
    if (!(s2 > 0)) stop("The estimated noise variance is not positive; reduce `r`.")
    W <- sweep(sv$v[, seq_len(r), drop = FALSE], 2, sqrt(pmax(kap[seq_len(r)] - s2, 0)), "*")
    unitvar <- s2
  } else {
    nz <- .pr_noise_structure(dat, "block", L, n, 8, NULL)
    est <- .pr_fit_ecme(t(X), r, nz, .pr_control(list()))
    W <- est$W
    if (ncol(W) < r) W <- cbind(W, matrix(0, D, r - ncol(W)))
    unitvar <- est$unitvar
  }
  psi <- unitvar[units$unit]
  Mw <- crossprod(W, W / psi)
  ev <- eigen((Mw + t(Mw)) / 2, symmetric = TRUE)
  W <- W %*% ev$vectors
  Wp <- W / psi
  M <- diag(r) + crossprod(W, Wp)
  Minv <- solve(M)
  list(W = W, E = Minv %*% t(Wp), lv = log(pmax(diag(Minv), 1e-12)), unitvar = unitvar,
       active = colSums(W^2) > 0)
}

# ------------------------------------------------------------------------------
# Parameters: flat named list of base-R arrays
# ------------------------------------------------------------------------------
.ae_spec <- function(D, r, arch, variational, units) {
  depth <- if (arch$width == 0) 0L else as.integer(arch$depth)
  list(D = D, r = r, width = as.integer(arch$width), depth = depth, activation = arch$activation,
       residual = isTRUE(arch$residual), variational = variational,
       nunits = if (variational) units$n else 0L, A = if (variational) units$A else NULL)
}

.ae_unif <- function(nr, nc, fan_in) {
  a <- 1 / sqrt(max(fan_in, 1))
  matrix(stats::runif(nr * nc, -a, a), nr, nc)
}

# lin: list(W, E, lv) or NULL (random linear part); seed: R RNG seed of the hidden layers
.ae_init_par <- function(spec, lin, seed, logpsi0 = NULL) {
  set.seed(seed)
  D <- spec$D
  r <- spec$r
  w <- spec$width
  par <- list()
  if (spec$residual) {
    if (is.null(lin)) {
      par$W <- .ae_unif(D, r, r)
      par$E <- .ae_unif(r, D, D)
    } else {
      par$W <- lin$W
      par$E <- lin$E
    }
  }
  if (spec$depth > 0) {
    fin <- r
    for (k in seq_len(spec$depth)) {
      par[[paste0("dec", k, "_w")]] <- .ae_unif(w, fin, fin)
      par[[paste0("dec", k, "_b")]] <- as.numeric(.ae_unif(w, 1, fin))
      fin <- w
    }
    par$dec_out <- if (spec$residual) matrix(0, D, w) else .ae_unif(D, w, w)
    fin <- D
    for (k in seq_len(spec$depth)) {
      par[[paste0("enc", k, "_w")]] <- .ae_unif(w, fin, fin)
      par[[paste0("enc", k, "_b")]] <- as.numeric(.ae_unif(w, 1, fin))
      fin <- w
    }
    par$enc_out <- if (spec$residual) matrix(0, r, w) else .ae_unif(r, w, w)
    if (spec$variational) par$enc_lv <- matrix(0, r, w)
  }
  if (spec$variational) {
    par$lv_b <- if (!is.null(lin) && !is.null(lin$lv)) as.numeric(lin$lv) else rep(0, r)
    par$logpsi <- if (!is.null(logpsi0)) as.numeric(logpsi0) else rep(0, spec$nunits)
  }
  par
}

.ae_is_linear_name <- function(nm) nm %in% c("W", "E")

# ------------------------------------------------------------------------------
# Torch forward maps and losses (functional; P = named list of tensors)
# ------------------------------------------------------------------------------
.ae_t_encode <- function(P, X, spec, act) {
  mu <- NULL
  if (!is.null(P$E)) mu <- torch::torch_matmul(X, P$E$t())
  lv <- NULL
  if (spec$depth > 0) {
    a <- X
    for (k in seq_len(spec$depth)) {
      a <- act(torch::torch_matmul(a, P[[paste0("enc", k, "_w")]]$t()) + P[[paste0("enc", k, "_b")]])
    }
    add <- torch::torch_matmul(a, P$enc_out$t())
    mu <- if (is.null(mu)) add else mu + add
    if (spec$variational) lv <- torch::torch_matmul(a, P$enc_lv$t()) + P$lv_b
  } else if (spec$variational) {
    lv <- P$lv_b$unsqueeze(1L)
  }
  list(mu = mu, lv = lv)
}

.ae_t_decode <- function(P, H, spec, act) {
  out <- NULL
  if (!is.null(P$W)) out <- torch::torch_matmul(H, P$W$t())
  if (spec$depth > 0) {
    a <- H
    for (k in seq_len(spec$depth)) {
      a <- act(torch::torch_matmul(a, P[[paste0("dec", k, "_w")]]$t()) + P[[paste0("dec", k, "_b")]])
    }
    add <- torch::torch_matmul(a, P$dec_out$t())
    out <- if (is.null(out)) add else out + add
  }
  out
}

# AESSA: mean squared error per coordinate; VAESSA: negative composite ELBO per
# window and coordinate (one reparametrised draw, or the supplied draws eps)
.ae_t_loss <- function(P, X, spec, act, Atens = NULL, eps = NULL) {
  q <- .ae_t_encode(P, X, spec, act)
  if (!spec$variational) {
    return(torch::torch_mean((X - .ae_t_decode(P, q$mu, spec, act))^2))
  }
  if (is.null(eps)) eps <- torch::torch_randn_like(q$mu)
  z <- q$mu + torch::torch_exp(0.5 * q$lv) * eps
  xh <- .ae_t_decode(P, z, spec, act)
  D <- spec$D
  lpsi <- torch::torch_matmul(Atens, P$logpsi)
  nll <- 0.5 * (D * log(2 * pi) + torch::torch_sum(lpsi) +
                  torch::torch_sum((X - xh)^2 / torch::torch_exp(lpsi), dim = 2L))
  kl <- 0.5 * torch::torch_sum(torch::torch_exp(q$lv) + q$mu^2 - 1 - q$lv, dim = 2L)
  torch::torch_mean(nll + kl) / D
}

.ae_dtype <- function(dtype) if (dtype == "float64") torch::torch_float64() else torch::torch_float32()

.ae_to_torch <- function(par, dtype) {
  lapply(par, function(p) torch::torch_tensor(p, dtype = .ae_dtype(dtype), requires_grad = TRUE))
}

.ae_from_torch <- function(P, template) {
  out <- lapply(names(template), function(nm) {
    v <- torch::as_array(P[[nm]]$detach()$to(dtype = torch::torch_float64()))
    if (is.matrix(template[[nm]])) matrix(v, nrow(template[[nm]]), ncol(template[[nm]])) else as.numeric(v)
  })
  names(out) <- names(template)
  out
}

# ------------------------------------------------------------------------------
# One training run (Adam, full batch or minibatches, early stopping on a
# validation block, best iterate kept). With monitor = FALSE the run lasts
# exactly `epochs` epochs and the last iterate is returned.
# ------------------------------------------------------------------------------
.ae_train_one <- function(Xtr, Xva, par0, spec, seed, epochs, lr, weight_decay, patience,
                          val_every, batch, dtype, floor_log = NULL, monitor = TRUE) {
  t0 <- proc.time()[["elapsed"]]
  set.seed(seed)
  torch::torch_manual_seed(seed)
  act <- .ae_act(spec$activation)$T
  dt <- .ae_dtype(dtype)
  P <- .ae_to_torch(par0, dtype)
  Atens <- if (spec$variational) torch::torch_tensor(spec$A, dtype = dt) else NULL
  Xt <- torch::torch_tensor(Xtr, dtype = dt)
  has_va <- !is.null(Xva) && nrow(Xva) > 0
  Xv <- if (has_va) torch::torch_tensor(Xva, dtype = dt) else NULL
  eps_tr <- eps_va <- NULL
  if (spec$variational) {
    eps_tr <- torch::torch_randn(nrow(Xtr), spec$r, dtype = dt)
    if (has_va) eps_va <- torch::torch_randn(nrow(Xva), spec$r, dtype = dt)
  }
  nms <- names(P)
  lin <- nms[.ae_is_linear_name(nms) | nms == "logpsi"]
  nl <- setdiff(nms, lin)
  groups <- list()
  if (length(lin)) groups[[length(groups) + 1]] <- list(params = unname(P[lin]), weight_decay = 0)
  if (length(nl)) groups[[length(groups) + 1]] <- list(params = unname(P[nl]), weight_decay = weight_decay)
  adam <- if (exists("optim_ignite_adam", where = asNamespace("torch"), inherits = FALSE)) {
    torch::optim_ignite_adam
  } else {
    torch::optim_adam
  }
  eval_loss <- function(X, eps) torch::with_no_grad(.ae_t_loss(P, X, spec, act, Atens, eps))$item()
  monitor_loss <- function() if (has_va) eval_loss(Xv, eps_va) else eval_loss(Xt, eps_tr)
  n <- nrow(Xtr)
  use_batch <- !is.null(batch) && batch < n
  tr_trace <- rep(NA_real_, epochs)
  va_ep <- integer(0)
  va_trace <- numeric(0)
  best <- if (monitor) monitor_loss() else NA_real_
  init_monitor <- best
  best_ep <- 0L
  best_state <- if (monitor) lapply(P, function(t) t$detach()$clone()) else NULL
  if (monitor) {
    va_ep <- 0L
    va_trace <- best
  }
  wait <- 0L
  ep <- 0L
  stopped <- "epochs"
  if (epochs > 0) {
    opt <- adam(groups, lr = lr)
    for (ep in seq_len(epochs)) {
      if (!use_batch) {
        opt$zero_grad()
        loss <- .ae_t_loss(P, Xt, spec, act, Atens)
        loss$backward()
        opt$step()
        tl <- loss$item()
      } else {
        perm <- sample.int(n)
        tl <- 0
        for (st in seq(1L, n, by = batch)) {
          ib <- perm[st:min(n, st + batch - 1L)]
          opt$zero_grad()
          loss <- .ae_t_loss(P, Xt[ib, , drop = FALSE], spec, act, Atens)
          loss$backward()
          opt$step()
          tl <- tl + loss$item() * length(ib) / n
        }
      }
      if (!is.null(floor_log)) {
        torch::with_no_grad(P$logpsi$copy_(torch::torch_maximum(P$logpsi,
                                                                torch::torch_tensor(floor_log, dtype = dt))))
      }
      tr_trace[ep] <- tl
      if (!is.finite(tl)) {
        stopped <- "non-finite loss"
        break
      }
      if (monitor && (ep %% val_every == 0L || ep == epochs)) {
        v <- monitor_loss()
        va_ep <- c(va_ep, ep)
        va_trace <- c(va_trace, v)
        if (!is.finite(v)) {
          stopped <- "non-finite validation loss"
          break
        }
        if (v < best - 1e-10 * (1 + abs(best))) {
          best <- v
          best_ep <- ep
          wait <- 0L
          best_state <- lapply(P, function(t) t$detach()$clone())
        } else {
          wait <- wait + val_every
          if (wait >= patience) {
            stopped <- "patience"
            break
          }
        }
      }
    }
  }
  if (monitor) {
    P <- best_state
  } else if (ep > 0 && !is.finite(tr_trace[ep])) {
    stop("training diverged (non-finite loss)")
  }
  list(par = .ae_from_torch(P, par0), best_epoch = if (monitor) best_ep else ep,
       best_monitor = best, init_monitor = init_monitor, epochs_run = ep, stopped = stopped,
       train_trace = tr_trace[seq_len(ep)],
       monitor_trace = data.frame(epoch = va_ep, loss = va_trace),
       time = proc.time()[["elapsed"]] - t0)
}

# ------------------------------------------------------------------------------
# Base-R forward maps on `par` (double precision)
# ------------------------------------------------------------------------------
.ae_mlp_R <- function(par, A, prefix, depth, act) {
  for (k in seq_len(depth)) {
    A <- act(sweep(A %*% t(par[[paste0(prefix, k, "_w")]]), 2, par[[paste0(prefix, k, "_b")]], "+"))
  }
  A
}

# encoder: X (n x D) -> list(mu = n x r, lv = n x r or NULL)
.ae_encode_R <- function(par, X, spec) {
  X <- as.matrix(X)
  n <- nrow(X)
  act <- .ae_act(spec$activation)$R
  mu <- if (!is.null(par$E)) X %*% t(par$E) else matrix(0, n, spec$r)
  lv <- NULL
  if (spec$depth > 0) {
    a <- .ae_mlp_R(par, X, "enc", spec$depth, act)
    mu <- mu + a %*% t(par$enc_out)
    if (spec$variational) lv <- sweep(a %*% t(par$enc_lv), 2, par$lv_b, "+")
  } else if (spec$variational) {
    lv <- matrix(par$lv_b, n, spec$r, byrow = TRUE)
  }
  list(mu = mu, lv = lv)
}

# decoder: H (n x r) -> n x D; part = "all", "linear" or "nonlinear"
.ae_decode_R <- function(par, H, spec, part = "all") {
  H <- matrix(H, ncol = spec$r)
  out <- matrix(0, nrow(H), spec$D)
  if (part != "nonlinear" && !is.null(par$W)) out <- out + H %*% t(par$W)
  if (part != "linear" && spec$depth > 0) {
    a <- .ae_mlp_R(par, H, "dec", spec$depth, .ae_act(spec$activation)$R)
    out <- out + a %*% t(par$dec_out)
  }
  out
}

# decoder value and Jacobian at one latent point h (length r): list(g = D, J = D x r)
.ae_jac_R <- function(par, h, spec) {
  r <- spec$r
  g <- numeric(spec$D)
  J <- matrix(0, spec$D, r)
  if (!is.null(par$W)) {
    g <- g + as.numeric(par$W %*% h)
    J <- J + par$W
  }
  if (spec$depth > 0) {
    ac <- .ae_act(spec$activation)
    a <- h
    M <- diag(r)
    for (k in seq_len(spec$depth)) {
      Wk <- par[[paste0("dec", k, "_w")]]
      pre <- as.numeric(Wk %*% a) + par[[paste0("dec", k, "_b")]]
      a <- ac$R(pre)
      M <- (ac$dR(pre, a) * Wk) %*% M
    }
    g <- g + as.numeric(par$dec_out %*% a)
    J <- J + par$dec_out %*% M
  }
  list(g = g, J = J)
}

# negative composite ELBO (per window, summed over windows) by Monte Carlo in base R
.ae_elbo_R <- function(par, X, spec, nsim = 20L, seed = 1L) {
  q <- .ae_encode_R(par, X, spec)
  n <- nrow(X)
  psi <- exp(as.numeric(spec$A %*% par$logpsi))
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed)
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
  ll <- 0
  for (s in seq_len(nsim)) {
    z <- q$mu + exp(0.5 * q$lv) * matrix(stats::rnorm(n * spec$r), n, spec$r)
    R <- X - .ae_decode_R(par, z, spec)
    ll <- ll + sum(-0.5 * (spec$D * log(2 * pi) + sum(log(psi)) + rowSums(sweep(R^2, 2, psi, "/"))))
  }
  kl_unit <- colMeans(0.5 * (exp(q$lv) + q$mu^2 - 1 - q$lv))
  list(elbo = ll / nsim - n * sum(kl_unit), kl_per_unit = kl_unit)
}
