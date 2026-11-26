#' @importFrom Rssa reconstruct
#' @export
Rssa::reconstruct

#' Reconstruction of an AESSA / VAESSA fit
#'
#' Diagonal averaging (Hankelization) of the decoded records
#' \eqn{g(h_j)}, \eqn{h_j = f(x_j)} (AESSA) or the posterior mean code (VAESSA),
#' split into the parts of the decoder: the linear part \eqn{W h_j} and the
#' nonlinear residual branch \eqn{\tilde g(h_j)}. The parts add up to the
#' signal. As in [prreconstruct()], the parts exclude the mean: the estimated
#' signal is `rec$mean + rec$signal` (use `+.funts` for functional input).
#'
#' @param x an object from [aessa()] or [vaessa()].
#' @param groups parts to return: any of `"signal"`, `"linear"`, `"nonlinear"`.
#' @param ... unused.
#' @return A list with the requested parts and `mean`, each in the class of the
#'   input (vector or `ts`, N x p matrix, list, or `funts`).
#' @seealso [aessa()], [prreconstruct()]
#' @examples
#' \donttest{
#' if (requireNamespace("torch", quietly = TRUE) && torch::torch_is_installed()) {
#'   set.seed(1)
#'   tt <- 1:150
#'   y <- tanh(3 * sin(2 * pi * tt / 25)) + rnorm(150, sd = 0.2)
#'   fit <- aessa(y, L = 50, r = 2, restarts = 1, epochs = 200, seed = 1)
#'   rec <- reconstruct(fit)
#'   plot(y, col = "grey")
#'   lines(rec$mean + rec$signal)
#' }
#' }
#' @method reconstruct aessa
#' @md
#' @export
reconstruct.aessa <- function(x, groups = c("signal", "linear", "nonlinear"), ...) {
  .ae_need_torch()
  groups <- match.arg(groups, c("signal", "linear", "nonlinear"), several.ok = TRUE)
  prep <- x$prep
  idx <- seq_len(prep$N)
  out <- list()
  for (g in groups) {
    part <- switch(g, signal = "all", linear = "linear", nonlinear = "nonlinear")
    Tm <- .ae_decode_R(x$par, x$scores, x$spec, part = part)
    out[[g]] <- .ae_records_to_input(prep, Tm, idx, add_mean = FALSE)
  }
  out$mean <- .ae_mean_input(prep)
  out
}

#' @rdname aeforecast
#' @param object an object from [aessa()] or [vaessa()].
#' @param ... arguments passed to [aeforecast()].
#' @importFrom stats predict
#' @md
#' @export
predict.aessa <- function(object, h = 1, ...) {
  aeforecast(object, h = h, ...)
}

.ae_model_name <- function(x) if (inherits(x, "vaessa")) "VAESSA" else "AESSA"

.ae_arch_text <- function(x) {
  a <- x$settings$arch
  if (a$width == 0) return("linear (no hidden layer)")
  sprintf("%s%s branch (%d hidden layer%s of width %d)",
          if (a$residual) "linear + " else "", a$activation, a$depth, if (a$depth > 1) "s" else "", a$width)
}

#' @export
print.aessa <- function(x, ...) {
  s <- x$settings
  cat(sprintf("%s fit (%s singular spectrum analysis)\n", .ae_model_name(x),
              if (inherits(x, "vaessa")) "variational autoencoder" else "autoencoder"))
  cat(sprintf("  N = %d, L = %d, K = %d windows, record dimension D = %d, latent dimension r = %d\n",
              s$N, s$L, s$K, s$D, s$r))
  cat(sprintf("  architecture: %s; start: %s\n", .ae_arch_text(x),
              if (s$init == "ssa") if (inherits(x, "vaessa")) "PrSSA" else "SSA" else "random"))
  cat(sprintf("  mean convention %s%s; %d restart(s), best = %d (epoch %d)%s\n", s$mean,
              if (isTRUE(s$scale)) ", scaled per block" else "", s$restarts, x$best_restart, x$best_epoch,
              if (!is.null(x$refit_trace)) ", refitted on all windows" else ""))
  b <- x$ssa_baseline
  lab <- if (inherits(x, "vaessa")) "negative ELBO per coordinate" else "MSE per coordinate"
  cat(sprintf("  %s (standardized scale): %.6g (linear start: %.6g, ratio %.4f)\n", lab, x$loss, b$loss,
              x$loss / b$loss))
  if (length(x$split$valid)) {
    cat(sprintf("  validation block (%d windows, gap %d): %.6g at start -> %.6g\n", length(x$split$valid),
                x$split$gap, b$valid_loss, b$final_valid_loss))
  }
  if (inherits(x, "vaessa")) {
    cat("  noise variance sigma2:", paste(format(x$sigma2, digits = 5), collapse = ", "),
        "  (PrSSA start:", paste(format(x$prssa_sigma2, digits = 5), collapse = ", "), ")\n")
    if (any(x$floor_active)) cat("  note: the sigma2 floor is attained\n")
    cat(sprintf("  active latent units (KL > 0.01): %d of %d\n", x$active_units, s$r))
  }
  cat(sprintf("  time: %.1f s\n", x$time))
  invisible(x)
}

#' @export
summary.aessa <- function(object, ...) {
  out <- list(model = .ae_model_name(object), settings = object$settings, restarts = object$restarts,
              best_restart = object$best_restart, best_epoch = object$best_epoch, loss = object$loss,
              ssa_baseline = object$ssa_baseline, sigma2 = object$sigma2,
              residual_ms = if (inherits(object, "vaessa")) object$residual_ms else object$sigma2,
              kl_per_unit = object$kl_per_unit, floor_active = object$floor_active,
              split = lengths(object$split[c("train", "valid")]), time = object$time)
  class(out) <- "summary.aessa"
  out
}

#' @export
print.summary.aessa <- function(x, ...) {
  s <- x$settings
  cat(sprintf("%s: N = %d, L = %d, r = %d, mean %s, lr = %g, weight decay = %g, dtype %s\n", x$model,
              s$N, s$L, s$r, s$mean, s$lr, s$weight_decay, s$dtype))
  cat(sprintf("windows: %d training, %d validation (gap %s)\n", x$split[1], x$split[2],
              if (x$split[2] > 0) s$gap else "-"))
  cat("restarts:\n")
  print(x$restarts, row.names = FALSE, digits = 6)
  cat(sprintf("final loss %.6g; linear start %.6g\n", x$loss, x$ssa_baseline$loss))
  cat("residual mean square per block (original scale):\n")
  print(x$residual_ms, digits = 6)
  if (!is.null(x$kl_per_unit)) {
    cat("noise variance:\n")
    print(x$sigma2, digits = 6)
    cat("KL per latent unit:", format(x$kl_per_unit, digits = 4), "\n")
  }
  invisible(x)
}

#' Plot an AESSA / VAESSA fit
#'
#' @param x an object from [aessa()] or [vaessa()].
#' @param type `"loss"` (training and monitored loss by epoch), `"latent"`
#'   (codes against the window index and, for r >= 2, the first two codes
#'   against each other) or `"components"` (data, signal, and its linear and
#'   nonlinear parts; first series for multivariate input).
#' @param ... further graphical parameters passed to `plot()`.
#' @return `x`, invisibly.
#' @md
#' @export
plot.aessa <- function(x, type = c("loss", "latent", "components"), ...) {
  type <- match.arg(type)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  if (type == "loss") {
    tr <- x$loss_trace
    va <- x$valid_trace
    rng <- range(c(tr, va$loss), finite = TRUE)
    plot(seq_along(tr), tr, type = "l", log = if (all(c(tr, va$loss) > 0, na.rm = TRUE)) "y" else "",
         ylim = rng, xlab = "epoch", ylab = "loss", main = paste(.ae_model_name(x), "training"), ...)
    if (nrow(va)) graphics::lines(va$epoch, va$loss, col = 2)
    graphics::abline(v = x$best_epoch, lty = 3)
    graphics::legend("topright", c("training", if (length(x$split$valid)) "validation" else "monitored"),
                     col = 1:2, lty = 1, bty = "n")
  } else if (type == "latent") {
    S <- x$scores
    if (ncol(S) >= 2) graphics::par(mfrow = c(1, 2))
    graphics::matplot(S, type = "l", lty = 1, xlab = "window j", ylab = "code", main = "latent codes", ...)
    if (ncol(S) >= 2) plot(S[, 1], S[, 2], type = "l", xlab = "h1", ylab = "h2", main = "h2 vs h1")
  } else {
    if (!x$prep$scalar) stop("type = \"components\" is available for vector and matrix input.")
    rec <- reconstruct(x)
    first <- function(v) if (is.matrix(v)) v[, 1] else as.numeric(v)
    y <- first(.pr_as_input(x$prep$dat, lapply(x$prep$dat$blocks, function(b) b$Z), seq_len(x$prep$N)))
    sig <- first(rec$signal) + first(rec$mean)
    graphics::par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))
    plot(y, type = "l", col = "grey50", ylab = "series", main = "data and signal", ...)
    graphics::lines(sig, col = 2)
    plot(first(rec$linear), type = "l", ylab = "linear part", main = "linear part W h")
    plot(first(rec$nonlinear), type = "l", ylab = "nonlinear part", main = "nonlinear branch")
  }
  invisible(x)
}

#' @export
print.aeforecast <- function(x, ...) {
  cat(sprintf("%s forecast, h = %d, method \"%s\"%s\n", x$model, x$h, x$method,
              if (x$interval != "none") sprintf(", %s intervals", x$interval) else ""))
  tab <- table(x$method_used)
  cat("  laws used:", paste(names(tab), tab, sep = " x ", collapse = ", "), "\n")
  if (!is.null(x$diagnostics)) {
    cat(sprintf("  decoder inversion: %d of %d horizons converged, max spread %.3g\n",
                sum(x$diagnostics$converged), nrow(x$diagnostics), max(x$diagnostics$spread, na.rm = TRUE)))
  }
  if (is.numeric(x$mean) && is.null(dim(x$mean))) {
    m <- cbind(mean = as.numeric(x$mean))
    if (!is.null(x$lower)) {
      for (nm in names(x$lower)) {
        m <- cbind(m, as.numeric(x$lower[[nm]]), as.numeric(x$upper[[nm]]))
        colnames(m)[ncol(m) - 1:0] <- paste(c("lo", "hi"), nm)
      }
    }
    print(utils::head(m, 12), digits = 5)
    if (nrow(m) > 12) cat("  ...\n")
  } else {
    cat("  point forecasts in $mean (", class(x$mean)[1], ")\n")
  }
  invisible(x)
}
