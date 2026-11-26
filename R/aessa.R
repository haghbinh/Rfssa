#' Network architecture of AESSA / VAESSA
#'
#' Describes the encoder and decoder networks used by [aessa()] and [vaessa()].
#' With `residual = TRUE` (the default) the decoder is \eqn{g(h) = W h +
#' \tilde g(h)} and the encoder is \eqn{f(x) = E x + \tilde f(x)}, where the
#' nonlinear branches \eqn{\tilde g, \tilde f} are multilayer perceptrons whose
#' output layers start at zero. The linear parts start at the Basic SSA
#' (AESSA) or PrSSA (VAESSA) solution, so the untrained network reproduces
#' SSA/PrSSA exactly and training can only move away from it when the loss
#' improves. With `width = 0` the model is linear.
#'
#' @param width number of units of each hidden layer of the nonlinear branches
#'   (0 gives a linear autoencoder; needs `residual = TRUE`).
#' @param depth number of hidden layers of each nonlinear branch.
#' @param activation activation function: `"tanh"` (default), `"sigmoid"` or `"softplus"`.
#' @param residual logical; if `TRUE`, linear part plus nonlinear residual branch
#'   (see Details). If `FALSE`, plain multilayer perceptrons with random
#'   initialization (the pilot architecture).
#' @return A list of class `"aessa_arch"`.
#' @seealso [aessa()], [vaessa()]
#' @examples
#' \donttest{
#' if (requireNamespace("torch", quietly = TRUE) && torch::torch_is_installed()) {
#'   aessa_arch(width = 16)
#' }
#' }
#' @md
#' @export
aessa_arch <- function(width = 32, depth = 1, activation = c("tanh", "sigmoid", "softplus"),
                       residual = TRUE) {
  .ae_need_torch()
  activation <- match.arg(activation)
  if (!is.numeric(width) || length(width) != 1 || width < 0 || width != round(width)) {
    stop("`width` must be a non-negative integer.")
  }
  if (!is.numeric(depth) || length(depth) != 1 || depth < 1 || depth != round(depth)) {
    stop("`depth` must be a positive integer.")
  }
  if (width == 0 && !isTRUE(residual)) stop("`width = 0` (linear model) needs `residual = TRUE`.")
  structure(list(width = as.integer(width), depth = as.integer(depth), activation = activation,
                 residual = isTRUE(residual)), class = "aessa_arch")
}

#' Autoencoder singular spectrum analysis (AESSA)
#'
#' Nonlinear SSA in which the rank-r SVD step of Basic SSA is replaced by an
#' autoencoder of the lag vectors (records). An encoder \eqn{f} maps each
#' record \eqn{x_j} (the L lagged values of every series or basis coordinate,
#' j = 1, ..., K) to an r-dimensional code \eqn{h_j}; a decoder \eqn{g} maps the
#' code back. The signal is the diagonal average (Hankelization) of the
#' decoded records, as in Basic SSA. With the default architecture
#' ([aessa_arch()], `residual = TRUE`) the decoder is linear part + nonlinear
#' residual, initialized at Basic SSA, so a fit that cannot improve on the
#' validation windows returns the SSA solution.
#'
#' @section Inputs: `Y` may be a numeric vector or `ts` (AESSA), a numeric matrix
#'   with one series per column (multivariate AESSA, records stacked series by
#'   series, as in PrMSSA/VMSSA), a [funts] object (functional AESSA on
#'   Gram-whitened basis coefficients, as in [prssa()]), or a list of these.
#'
#' @section Training: The loss is the mean squared error of the records (the
#'   window loss; its linear minimizer is Basic SSA). The windows are split into
#'   a training block (the first windows) and a validation block (the last
#'   `valid` fraction), separated by `gap` windows so that no observation is
#'   shared. Each of the `restarts` runs starts from the linear solution of the
#'   training windows with a different random initialization of the hidden
#'   layers, is trained by Adam (`torch::optim_ignite_adam` when available) with
#'   weight decay on the nonlinear branches only, and keeps the iterate with the
#'   smallest validation loss (epoch 0 = the linear solution). If
#'   `refit = TRUE`, the best restart is then re-run on all windows, from the
#'   linear solution of all windows and the same seed, for its selected number
#'   of epochs. With `valid = 0`, the runs use all windows and early stopping
#'   monitors the training loss.
#'
#' @param Y the data: see Section Inputs.
#' @param L window length, 2 <= L <= N - 1.
#' @param r latent dimension (bottleneck), 1 <= r < min(D, K).
#' @param mean mean convention: `"M1"` (default) centres each series at its
#'   sample mean, `"M0"` uses no centring (Basic SSA), `"M2"` centres the
#'   records at their mean record (row-centred SSA).
#' @param scale logical; divide each block (series or functional variable) by its
#'   root mean square before fitting (undone in all outputs).
#' @param arch network architecture, from [aessa_arch()].
#' @param init `"ssa"` (default): linear part initialized at SSA (AESSA) or
#'   PrSSA (VAESSA); `"random"`: random linear part.
#' @param weight_decay weight decay (L2 penalty in Adam) of the nonlinear branches.
#' @param lr Adam learning rate.
#' @param epochs maximum number of epochs per run (0 returns the initial network).
#' @param batch minibatch size (`NULL`: full batch).
#' @param patience early stopping patience, in epochs.
#' @param valid fraction of windows in the validation block (0: none).
#' @param gap number of windows between the training and validation blocks
#'   (default `L - 1`: no shared observation).
#' @param refit logical; re-run the best restart on all windows (see Section Training).
#' @param restarts number of restarts (random initializations of the hidden layers).
#' @param seed integer seed (`NULL`: drawn from the R random number generator).
#' @param dtype `"float64"` (default) or `"float32"` torch computations.
#' @param val_every evaluate the monitored loss every `val_every` epochs.
#' @param verbose logical; print one line per restart.
#' @return An object of class `"aessa"`: a list with
#'   \item{par}{network parameters (named list of matrices and vectors; `W`, `E` are the linear parts)}
#'   \item{scores}{K x r matrix of codes \eqn{h_j = f(x_j)}}
#'   \item{records_hat}{K x D matrix of decoded records (standardized, centred working scale)}
#'   \item{signal}{the reconstructed signal in the class of `Y` (mean included)}
#'   \item{loss}{training-window loss of the final network (mean squared error per coordinate, standardized scale)}
#'   \item{ssa_baseline}{list with `loss` (the same loss of the linear start on all windows), `valid_loss` and `final_valid_loss` (validation loss at epoch 0 and at the selected epoch of the best restart)}
#'   \item{restarts}{data frame with one row per restart (seed, best epoch, epochs run, losses, time)}
#'   \item{best_restart, best_epoch}{the selected restart and epoch}
#'   \item{loss_trace, valid_trace}{training loss per epoch and monitored loss of the selected restart (stage 1)}
#'   \item{refit_trace}{training loss per epoch of the refit (or `NULL`)}
#'   \item{sigma2}{residual mean square of the records per block (original scale)}
#'   \item{prep, spec, settings}{data layout, network specification and the arguments}
#' @references Hsieh, W. W. and Wu, A. (2002). Nonlinear multichannel singular
#'   spectrum analysis of the tropical Pacific climate variability using a
#'   neural network approach. \emph{Journal of Geophysical Research}, 107(C7), 3076.
#' @seealso [vaessa()], [aeforecast()], [reconstruct.aessa()], [aessa_select()], [prssa()]
#' @examples
#' \donttest{
#' if (requireNamespace("torch", quietly = TRUE) && torch::torch_is_installed()) {
#'   set.seed(1)
#'   tt <- 1:200
#'   y <- tanh(3 * sin(2 * pi * tt / 25)) + rnorm(200, sd = 0.2)
#'   fit <- aessa(y, L = 50, r = 2, restarts = 1, epochs = 300, seed = 1)
#'   fit
#'   fc <- aeforecast(fit, h = 10)
#'   fc$mean
#' }
#' }
#' @md
#' @export
aessa <- function(Y, L, r, mean = c("M1", "M0", "M2"), scale = TRUE, arch = aessa_arch(),
                  init = c("ssa", "random"), weight_decay = 1e-3, lr = 1e-3, epochs = 3000,
                  batch = NULL, patience = 200, valid = 0.15, gap = NULL, refit = TRUE,
                  restarts = 3, seed = NULL, dtype = c("float64", "float32"), val_every = 5,
                  verbose = FALSE) {
  .ae_need_torch()
  .ae_fit(Y = Y, L = L, r = r, variational = FALSE, noise = "iso", mean = match.arg(mean),
          scale = scale, arch = arch, init = match.arg(init), weight_decay = weight_decay, lr = lr,
          epochs = epochs, batch = batch, patience = patience, valid = valid, gap = gap,
          refit = refit, restarts = restarts, seed = seed, dtype = match.arg(dtype),
          sigma_floor = NA_real_, val_every = val_every, verbose = verbose, call = match.call())
}

#' Variational autoencoder singular spectrum analysis (VAESSA)
#'
#' The probabilistic (variational) counterpart of [aessa()] and the nonlinear
#' counterpart of PrSSA ([prssa()]). Every record is modelled as
#' \eqn{x_j = g(h_j) + \epsilon_j}, \eqn{h_j \sim N(0, I_r)}, \eqn{\epsilon_j
#' \sim N(0, \Psi)}, with \eqn{\Psi} isotropic (`noise = "iso"`) or constant
#' within each block (`noise = "block"`, one variance per series or functional
#' variable). The windows overlap, so the sum of the per-window evidence lower
#' bounds is a composite ELBO, maximized over the decoder, \eqn{\Psi} and an
#' amortized diagonal-Gaussian encoder \eqn{q(h_j | x_j)} by stochastic
#' gradients (one reparametrized draw per window and epoch). With a linear
#' decoder the optimum is PrSSA; the default architecture starts exactly at the
#' PrSSA fit (loadings, noise variances and the exact posterior as encoder).
#'
#' Each noise variance is bounded below by `sigma_floor` times its PrSSA value
#' (a safeguard against the unbounded likelihood of flexible decoders);
#' `floor_active` reports whether the bound is attained. Training, validation
#' and refitting are as in [aessa()]; the monitored loss is the negative ELBO
#' with fixed draws.
#'
#' @inheritParams aessa
#' @param noise `"iso"` (default) or `"block"`.
#' @param sigma_floor lower bound for each noise variance, as a fraction of its
#'   PrSSA estimate (0 disables it).
#' @return An object of class `c("vaessa", "aessa")`: the elements of [aessa()]
#'   (with `scores` = posterior means of the codes and `loss` = negative ELBO per
#'   window and coordinate) plus
#'   \item{sigma2}{noise variance(s) on the original scale (a number for `"iso"`, a named vector for `"block"`)}
#'   \item{psi}{diagonal of \eqn{\Psi} on the standardized scale (length D)}
#'   \item{elbo}{composite ELBO of all windows on the scale of the data (Monte Carlo, 20 draws per window); comparable with the composite log-likelihood `prssa()$loglik`}
#'   \item{kl_per_unit}{average KL divergence per latent unit; units with KL < 0.01 nats are inactive}
#'   \item{prssa_sigma2}{PrSSA noise variance(s) of the linear start (original scale)}
#'   \item{floor_active}{logical, per noise unit}
#'   \item{scores_var}{K x r posterior variances of the codes}
#' @seealso [aessa()], [aeforecast()], [prssa()]
#' @examples
#' \donttest{
#' if (requireNamespace("torch", quietly = TRUE) && torch::torch_is_installed()) {
#'   set.seed(1)
#'   tt <- 1:200
#'   y <- tanh(3 * sin(2 * pi * tt / 25)) + rnorm(200, sd = 0.2)
#'   fit <- vaessa(y, L = 50, r = 2, restarts = 1, epochs = 300, seed = 1)
#'   fit$sigma2
#'   fc <- aeforecast(fit, h = 10, interval = "laplace")
#'   cbind(fc$mean, fc$lower[[1]], fc$upper[[1]])
#' }
#' }
#' @md
#' @export
vaessa <- function(Y, L, r, noise = c("iso", "block"), mean = c("M1", "M0", "M2"), scale = TRUE,
                   arch = aessa_arch(), init = c("ssa", "random"), weight_decay = 1e-3, lr = 1e-3,
                   epochs = 3000, batch = NULL, patience = 200, valid = 0.15, gap = NULL,
                   refit = TRUE, restarts = 3, seed = NULL, dtype = c("float64", "float32"),
                   sigma_floor = 0.25, val_every = 5, verbose = FALSE) {
  .ae_need_torch()
  .ae_fit(Y = Y, L = L, r = r, variational = TRUE, noise = match.arg(noise), mean = match.arg(mean),
          scale = scale, arch = arch, init = match.arg(init), weight_decay = weight_decay, lr = lr,
          epochs = epochs, batch = batch, patience = patience, valid = valid, gap = gap,
          refit = refit, restarts = restarts, seed = seed, dtype = match.arg(dtype),
          sigma_floor = sigma_floor, val_every = val_every, verbose = verbose, call = match.call())
}

# ------------------------------------------------------------------------------
# Fitting core shared by aessa() and vaessa()
# ------------------------------------------------------------------------------
.ae_fit <- function(Y, L, r, variational, noise, mean, scale, arch, init, weight_decay, lr, epochs,
                    batch, patience, valid, gap, refit, restarts, seed, dtype, sigma_floor,
                    val_every, verbose, call) {
  t0 <- proc.time()[["elapsed"]]
  if (!inherits(arch, "aessa_arch")) stop("`arch` must be created by aessa_arch().")
  if (init == "ssa" && !arch$residual) stop("`init = \"ssa\"` needs `arch$residual = TRUE`.")
  chk_num <- function(x, nm, lo = 0) {
    if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < lo) {
      stop(sprintf("`%s` must be a number >= %s.", nm, lo))
    }
  }
  chk_num(weight_decay, "weight_decay")
  chk_num(lr, "lr", 1e-12)
  chk_num(epochs, "epochs")
  chk_num(patience, "patience", 1)
  chk_num(restarts, "restarts", 1)
  chk_num(val_every, "val_every", 1)
  if (!is.numeric(valid) || length(valid) != 1 || valid < 0 || valid >= 0.5) {
    stop("`valid` must be a number in [0, 0.5).")
  }
  epochs <- as.integer(epochs)
  restarts <- as.integer(restarts)
  val_every <- as.integer(val_every)
  if (!is.null(batch)) chk_num(batch, "batch", 1)
  prep <- .ae_prepare(Y, L, mean, scale)
  K <- prep$K
  D <- prep$D
  L <- prep$L
  if (missing(r) || length(r) != 1 || !is.numeric(r) || r != round(r) || r < 1 || r >= min(D, K)) {
    stop("`r` must be an integer with 1 <= r < min(D, K).")
  }
  r <- as.integer(r)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max %/% 4L, 1L)
  if (!is.numeric(seed) || length(seed) != 1) stop("`seed` must be a single integer.")
  seed <- as.integer(seed)
  units <- .ae_units(prep, noise)
  spec <- .ae_spec(D, r, arch, variational, units)
  X <- prep$X
  # ---- training / validation split -----------------------------------------
  if (is.null(gap)) gap <- L - 1L
  gap <- as.integer(gap)
  nmin <- max(10L, 2L * r + 2L)
  if (valid > 0) {
    nv <- max(1L, as.integer(round(valid * K)))
    if (K - nv - gap < nmin) {
      gap_new <- max(0L, K - nv - nmin)
      if (K - nv - gap_new < nmin) stop("Too few windows for a validation block: use `valid = 0` or a smaller `L`.")
      warning(sprintf("Gap between training and validation windows reduced from %d to %d (few windows).",
                      gap, gap_new))
      gap <- gap_new
    }
    tr_idx <- seq_len(K - nv - gap)
    va_idx <- seq.int(K - nv + 1L, K)
  } else {
    tr_idx <- seq_len(K)
    va_idx <- integer(0)
  }
  # ---- linear starting values ------------------------------------------------
  lin_start <- function(Xs) {
    if (!variational) {
      l <- .ae_linear_ssa(Xs, r)
      return(list(W = l$W, E = l$E))
    }
    .ae_linear_ppca(Xs, r, units, prep$dat, L, noise)
  }
  lin1 <- lin_start(X[tr_idx, , drop = FALSE])
  lin_all <- if (length(va_idx)) lin_start(X) else lin1
  floor_log <- NULL
  logpsi1 <- logpsi_all <- NULL
  if (variational) {
    logpsi1 <- log(lin1$unitvar)
    logpsi_all <- log(lin_all$unitvar)
    if (is.numeric(sigma_floor) && length(sigma_floor) == 1 && sigma_floor > 0) {
      floor_log <- log(sigma_floor) + logpsi_all
    }
  }
  lin_use1 <- if (init == "ssa") lin1 else if (variational) list(lv = NULL) else NULL
  if (init == "random") lin_use1 <- NULL
  # ---- stage 1: restarts ----------------------------------------------------
  Xtr <- X[tr_idx, , drop = FALSE]
  Xva <- if (length(va_idx)) X[va_idx, , drop = FALSE] else NULL
  runs <- vector("list", restarts)
  seeds <- seed + seq_len(restarts) - 1L
  errs <- character(0)
  for (i in seq_len(restarts)) {
    par0 <- .ae_init_par(spec, lin_use1, seeds[i], logpsi0 = logpsi1)
    runs[[i]] <- tryCatch(
      .ae_train_one(Xtr, Xva, par0, spec, seeds[i], epochs, lr, weight_decay, patience, val_every,
                    batch, dtype, floor_log = floor_log, monitor = TRUE),
      error = function(e) {
        errs <<- c(errs, conditionMessage(e))
        NULL
      })
    if (verbose && !is.null(runs[[i]])) {
      message(sprintf("restart %d: best epoch %d of %d, monitored loss %.6g (start %.6g), %.1f s",
                      i, runs[[i]]$best_epoch, runs[[i]]$epochs_run, runs[[i]]$best_monitor,
                      runs[[i]]$init_monitor, runs[[i]]$time))
    }
  }
  ok <- !vapply(runs, is.null, logical(1))
  if (!any(ok)) stop("All restarts failed: ", paste(unique(errs), collapse = "; "))
  mon <- vapply(runs, function(u) if (is.null(u)) Inf else u$best_monitor, numeric(1))
  b <- which.min(mon)
  best <- runs[[b]]
  # ---- stage 2: refit on all windows ------------------------------------------
  refit_trace <- NULL
  if (length(va_idx) && isTRUE(refit)) {
    lin_use2 <- if (init == "ssa") lin_all else NULL
    par0 <- .ae_init_par(spec, lin_use2, seeds[b], logpsi0 = logpsi_all)
    if (best$best_epoch > 0) {
      rf <- .ae_train_one(X, NULL, par0, spec, seeds[b], best$best_epoch, lr, weight_decay, patience,
                          val_every, batch, dtype, floor_log = floor_log, monitor = FALSE)
      par <- rf$par
      refit_trace <- rf$train_trace
    } else {
      par <- par0
      refit_trace <- numeric(0)
    }
  } else {
    par <- best$par
  }
  # ---- outputs ------------------------------------------------------------------
  enc <- .ae_encode_R(par, X, spec)
  Xh <- .ae_decode_R(par, enc$mu, spec)
  signal <- .ae_records_to_input(prep, Xh, seq_len(prep$N))
  resid <- X - Xh
  rms_block <- vapply(seq_along(prep$dims), function(bb) {
    mean(resid[, prep$blk == bb]^2) * prep$scale[bb]^2
  }, numeric(1))
  names(rms_block) <- prep$dat$names
  # the linear (SSA / PrSSA) fit of all windows, for the baseline losses
  spec_lin <- spec
  spec_lin$depth <- 0L
  spec_lin$width <- 0L
  spec_lin$residual <- TRUE
  lin_par <- .ae_init_par(spec_lin, lin_all, seed, logpsi0 = logpsi_all)
  lin_enc <- .ae_encode_R(lin_par, X, spec_lin)
  restarts_df <- data.frame(
    restart = seq_len(restarts), seed = seeds,
    best_epoch = vapply(runs, function(u) if (is.null(u)) NA_real_ else u$best_epoch, numeric(1)),
    epochs_run = vapply(runs, function(u) if (is.null(u)) NA_real_ else u$epochs_run, numeric(1)),
    start_loss = vapply(runs, function(u) if (is.null(u)) NA_real_ else u$init_monitor, numeric(1)),
    best_loss = ifelse(is.finite(mon), mon, NA_real_),
    stopped = vapply(runs, function(u) if (is.null(u)) "failed" else u$stopped, character(1)),
    time = vapply(runs, function(u) if (is.null(u)) NA_real_ else u$time, numeric(1))
  )
  fit <- list(par = par, scores = enc$mu, records_hat = Xh, signal = signal,
              restarts = restarts_df, best_restart = b, best_epoch = best$best_epoch,
              loss_trace = best$train_trace, valid_trace = best$monitor_trace,
              refit_trace = refit_trace, errors = errs,
              split = list(train = tr_idx, valid = va_idx, gap = if (length(va_idx)) gap else NA_integer_),
              prep = prep, spec = spec, units = units,
              settings = list(L = L, r = r, N = prep$N, K = K, D = D, mean = mean, scale = scale,
                              arch = arch, init = init, weight_decay = weight_decay, lr = lr,
                              epochs = epochs, batch = batch, patience = patience, valid = valid,
                              gap = gap, refit = refit, restarts = restarts, seed = seed,
                              dtype = dtype, val_every = val_every, noise = noise,
                              sigma_floor = sigma_floor, variational = variational),
              call = call)
  if (!variational) {
    fit$loss <- mean(resid^2)
    fit$ssa_baseline <- list(loss = mean((X - .ae_decode_R(lin_par, lin_enc$mu, spec_lin))^2),
                             valid_loss = best$init_monitor, final_valid_loss = best$best_monitor)
    fit$sigma2 <- rms_block
    class(fit) <- "aessa"
  } else {
    psi <- exp(as.numeric(units$A %*% par$logpsi))
    el <- .ae_elbo_R(par, X, spec, nsim = 20L, seed = seed)
    el0 <- .ae_elbo_R(lin_par, X, spec_lin, nsim = 20L, seed = seed)
    s2u <- exp(par$logpsi)
    to_orig <- function(v) {
      if (noise == "iso" && length(unique(prep$scale)) == 1) return(v * prep$scale[1]^2)
      out <- if (noise == "iso") v * prep$scale^2 else v * prep$scale^2
      names(out) <- prep$dat$names
      out
    }
    s2u_orig <- to_orig(s2u)
    pr_orig <- to_orig(exp(logpsi_all))
    # ELBO on the scale of the data (Jacobian of the standardization), comparable with prssa()$loglik
    logjac <- K * sum(log(prep$scale[prep$blk]))
    fit$loss <- -el$elbo / (K * D)
    fit$elbo <- el$elbo - logjac
    fit$kl_per_unit <- el$kl_per_unit
    fit$active_units <- sum(el$kl_per_unit > 0.01)
    fit$sigma2 <- s2u_orig
    fit$sigma2_std <- s2u
    fit$psi <- psi
    fit$prssa_sigma2 <- pr_orig
    fit$floor_active <- if (is.null(floor_log)) rep(FALSE, units$n) else abs(par$logpsi - floor_log) < 1e-8
    fit$scores_var <- exp(enc$lv)
    fit$residual_ms <- rms_block
    fit$ssa_baseline <- list(loss = -el0$elbo / (K * D), elbo = el0$elbo - logjac,
                             valid_loss = best$init_monitor, final_valid_loss = best$best_monitor)
    class(fit) <- c("vaessa", "aessa")
  }
  fit$time <- proc.time()[["elapsed"]] - t0
  fit
}
