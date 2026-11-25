# =============================================================================
# S3 methods for "prssa", "prforecast" and "prrank" objects.
# =============================================================================

.pr_model_name <- function(x) {
  dat <- x$dat
  kinds <- vapply(dat$blocks, function(b) b$kind, character(1))
  base <- if (all(kinds == "scalar")) {
    if (length(kinds) == 1) "PrSSA" else "PrMSSA"
  } else if (all(kinds == "functional")) {
    "PrFSSA"
  } else {
    "PrHSSA (hybrid)"
  }
  if (x$alpha > 0) base <- paste0(base, " (smoothed, alpha = ", format(x$alpha), ")")
  base
}

#' Print, Summary, Log-Likelihood and Plot Methods for PrSSA Fits
#'
#' Methods for objects of class \code{"prssa"} returned by \code{\link{prssa}}.
#'
#' @param x,object An object of class \code{"prssa"} (or \code{"summary.prssa"}).
#' @param digits Number of significant digits printed.
#' @param n Number of eigenvalues shown in the scree plot.
#' @param ... Further arguments (passed to \code{plot} for the plot method).
#'
#' @return \code{print} returns its argument invisibly. \code{summary} returns an object of
#'   class \code{"summary.prssa"} with a component table (eigenvalues, shrinkage), the noise
#'   estimates, the composite log-likelihood, the nominal number of parameters, the nominal
#'   CL-BIC and, for isotropic noise without smoothing, an approximate 95\% interval for
#'   \eqn{\sigma^2} based on an overlap-aware HAC variance. \code{logLik} returns a
#'   \code{"logLik"} object with attributes \code{df} (nominal parameter count) and
#'   \code{nobs} (number of windows \eqn{K}). \code{plot} draws a scree plot of the leading
#'   eigenvalues on a log scale with the noise level as a horizontal line and returns
#'   \code{NULL} invisibly.
#'
#' @examples
#' set.seed(1)
#' tt <- 1:200
#' y <- 2 * sin(2 * pi * tt / 20) + cos(2 * pi * tt / 7) + rnorm(200, sd = 0.7)
#' fit <- prssa(y, L = 24, r = 4)
#' summary(fit)
#' logLik(fit)
#' plot(fit)
#'
#' @name prssa-methods
NULL

#' @rdname prssa-methods
#' @export
print.prssa <- function(x, digits = 4, ...) {
  cat(.pr_model_name(x), "fit:", x$noise, "noise,", x$mean_type, "mean\n")
  cat("  N =", x$N, " L =", x$L, " K =", x$K, " record dimension D =", x$D,
      " rank r =", x$r, "\n")
  s2 <- x$sigma2
  if (is.list(s2)) {
    cat("  noise variance (range per block):\n")
    for (nm in names(s2)) cat("   ", nm, ":", paste(format(range(s2[[nm]]), digits = digits), collapse = " - "), "\n")
  } else {
    cat("  noise variance:", paste(names(s2), format(s2, digits = digits), sep = if (is.null(names(s2))) "" else " = ", collapse = ", "), "\n")
  }
  cat("  shrinkage omega:", format(x$omega, digits = digits), "\n")
  cat("  composite log-likelihood:", format(x$loglik, digits = digits + 4),
      " (nominal df =", format(x$df), ")\n")
  if (x$noise != "iso") cat("  ECME iterations:", x$iterations, if (x$converged) "(converged)" else "(NOT converged)", "\n")
  invisible(x)
}

#' @rdname prssa-methods
#' @export
summary.prssa <- function(object, ...) {
  tab <- data.frame(component = seq_len(object$r),
                    eigenvalue = object$eigenvalues[seq_len(object$r)],
                    lambda = object$lambda, omega = object$omega)
  ci <- NULL
  if (object$noise == "iso" && object$alpha == 0) ci <- .pr_sigma2_ci(object)
  out <- list(model = .pr_model_name(object), noise = object$noise, mean = object$mean_type,
              N = object$N, L = object$L, K = object$K, D = object$D, r = object$r,
              table = tab, sigma2 = object$sigma2, sigma2_ci = ci,
              loglik = object$loglik, df = object$df,
              clbic_nominal = -2 * object$loglik + object$df * log(object$K),
              iterations = object$iterations, converged = object$converged)
  class(out) <- "summary.prssa"
  out
}

#' @rdname prssa-methods
#' @export
print.summary.prssa <- function(x, digits = 4, ...) {
  cat(x$model, "with", x$noise, "noise and", x$mean, "mean\n")
  cat("N =", x$N, " L =", x$L, " K =", x$K, " D =", x$D, " r =", x$r, "\n\n")
  print(format(x$table, digits = digits), row.names = FALSE)
  cat("\n")
  if (is.list(x$sigma2)) {
    cat("noise variance (range per block):\n")
    for (nm in names(x$sigma2)) cat("  ", nm, ":", paste(format(range(x$sigma2[[nm]]), digits = digits), collapse = " - "), "\n")
  } else {
    cat("noise variance:", format(x$sigma2, digits = digits), "\n")
  }
  if (!is.null(x$sigma2_ci)) {
    cat("approximate 95% interval for sigma2 (overlap-aware HAC):",
        paste(format(x$sigma2_ci, digits = digits), collapse = " - "), "\n")
  }
  cat("composite log-likelihood:", format(x$loglik, digits = digits + 4), "\n")
  cat("nominal df:", x$df, "  nominal CL-BIC:", format(x$clbic_nominal, digits = digits + 4), "\n")
  invisible(x)
}

#' @rdname prssa-methods
#' @importFrom stats logLik
#' @export
logLik.prssa <- function(object, ...) {
  structure(object$loglik, df = object$df, nobs = object$K, class = "logLik")
}

#' @rdname prssa-methods
#' @export
plot.prssa <- function(x, n = NULL, ...) {
  ev <- x$eigenvalues
  if (is.null(n)) n <- length(ev)
  n <- min(n, length(ev))
  ylab <- if (x$noise == "iso" && x$alpha == 0) "eigenvalue" else "noise-whitened eigenvalue"
  graphics::plot(seq_len(n), ev[seq_len(n)], log = "y", type = "b", pch = 19,
                 xlab = "component", ylab = ylab, main = .pr_model_name(x), ...)
  graphics::abline(h = x$eig_scale, lty = 2, col = "red")
  graphics::abline(v = x$r + 0.5, lty = 3, col = "grey40")
  invisible(NULL)
}

# approximate interval for sigma2 (isotropic noise): mean of the per-window
# residual energies with an overlap-aware HAC long-run variance
.pr_sigma2_ci <- function(fit, level = 0.95) {
  X <- .pr_fit_trajectory(fit, smoothed = FALSE)
  K <- ncol(X)
  D <- nrow(X)
  z <- (colSums(X^2) - colSums(crossprod(fit$U, X)^2)) / (D - fit$r)
  zc <- z - mean(z)
  M0 <- floor(4 * (K / 100)^(2 / 9))
  gam <- function(m) sum(zc[1:(K - m)] * zc[(1 + m):K]) / K
  mx <- min(fit$L - 1 + M0, K - 1)
  g <- vapply(seq_len(mx), gam, numeric(1))
  w_ov <- c(rep(1, min(fit$L - 1, mx)), 1 - seq_len(max(0, mx - (fit$L - 1))) / (M0 + 1))
  w_ba <- 1 - seq_len(mx) / (mx + 1)
  lrv <- sum(zc^2) / K + 2 * sum(w_ov[seq_len(mx)] * g)
  if (!(lrv > 0)) lrv <- sum(zc^2) / K + 2 * sum(w_ba * g)
  zq <- stats::qnorm(1 - (1 - level) / 2)
  mean(z) + c(-1, 1) * zq * sqrt(max(lrv, 0) / K)
}

#' @rdname prssa-methods
#' @export
print.prforecast <- function(x, ...) {
  cat("PrSSA forecast:", x$h, "step(s), method", x$method)
  if (any(x$method_used != x$method)) cat(" (", paste(unique(x$method_used), collapse = " + "), ")", sep = "")
  cat("\n")
  if (is.numeric(x$mean) && is.null(dim(x$mean))) {
    tab <- data.frame(h = seq_len(x$h), mean = as.numeric(x$mean), sd = x$sd)
    for (nm in names(x$lower)) {
      tab[[paste("lower", nm)]] <- as.numeric(x$lower[[nm]])
      tab[[paste("upper", nm)]] <- as.numeric(x$upper[[nm]])
    }
    print(tab, row.names = FALSE, digits = 4)
  } else {
    cat("levels:", paste(format(x$level), collapse = ", "), " bands:", paste(x$bands, collapse = ", "), "\n")
    cat("components: mean, sd, lower, upper, intervals, cov\n")
  }
  invisible(x)
}

#' @rdname prssa-methods
#' @export
print.prrank <- function(x, ...) {
  cat("PrSSA rank selection; criterion:", x$criterion, "; selected r =", x$selected, "\n")
  print(x$table, row.names = FALSE, digits = 6)
  invisible(x)
}
