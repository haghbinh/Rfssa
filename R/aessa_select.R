#' Latent dimension selection for AESSA / VAESSA by rolling-origin validation
#'
#' For each candidate latent dimension, fits [aessa()] or [vaessa()] to the
#' data up to each of `origins` forecast origins (the last one at N - h), forecasts
#' `h` steps with [aeforecast()] and records the mean squared forecast error.
#' The selected r minimizes the error averaged over origins and horizons.
#' Rolling-origin validation is more stable than a single split (see the AESSA
#' pilot study); the PrSSA criterion [prrank()] gives the rank of the linear
#' model, a natural upper end for `r`.
#'
#' @param Y numeric vector or `ts`, or numeric matrix (one series per column).
#' @param L window length.
#' @param r candidate latent dimensions.
#' @param model `"aessa"` or `"vaessa"`.
#' @param h forecast horizon of the validation (default: the smaller of `L - max(r)` and half of `L`, at least 1).
#' @param origins number of forecast origins.
#' @param step distance between consecutive origins (default `h`).
#' @param method forecasting method passed to [aeforecast()].
#' @param ... further arguments of [aessa()] / [vaessa()] (for example `restarts`, `seed`, `epochs`).
#' @param verbose logical; report progress.
#' @return An object of class `"aessa_select"`: `r` (selected), `table` (mean
#'   squared error by r), `details` (one row per r and origin), and the settings.
#' @seealso [aessa()], [vaessa()], [prrank()]
#' @examples
#' \donttest{
#' if (requireNamespace("torch", quietly = TRUE) && torch::torch_is_installed()) {
#'   set.seed(3)
#'   tt <- 1:160
#'   y <- tanh(3 * sin(2 * pi * tt / 20)) + rnorm(160, sd = 0.2)
#'   sel <- aessa_select(y, L = 40, r = 1:2, origins = 2, h = 10, restarts = 1,
#'                       epochs = 100, seed = 1)
#'   sel
#' }
#' }
#' @md
#' @export
aessa_select <- function(Y, L, r = 1:4, model = c("aessa", "vaessa"), h = NULL, origins = 3,
                         step = NULL, method = c("window", "recursive", "latent"), ...,
                         verbose = FALSE) {
  .ae_need_torch()
  model <- match.arg(model)
  method <- match.arg(method)
  is_mat <- is.matrix(Y) || is.data.frame(Y)
  if (is_mat) {
    Y <- as.matrix(Y)
  } else if (is.numeric(Y)) {
    Y <- as.numeric(Y)
  } else {
    stop("`Y` must be a numeric vector, ts, or matrix.")
  }
  N <- if (is_mat) nrow(Y) else length(Y)
  if (!is.numeric(r) || any(r < 1) || any(r != round(r))) stop("`r` must contain positive integers.")
  r <- sort(unique(as.integer(r)))
  if (is.null(h)) h <- min(L - max(r), max(1L, L %/% 2L))
  if (h < 1) stop("`h` must be at least 1.")
  h <- as.integer(h)
  if (is.null(step)) step <- h
  ends <- N - h - (rev(seq_len(origins)) - 1L) * step
  if (min(ends) - L + 1L < 20L) stop("Too few observations for these origins: reduce `origins`, `step` or `L`.")
  fitfun <- if (model == "aessa") aessa else vaessa
  rows <- list()
  for (rr in r) {
    for (o in seq_along(ends)) {
      n0 <- ends[o]
      Yfit <- if (is_mat) Y[seq_len(n0), , drop = FALSE] else Y[seq_len(n0)]
      Ytest <- if (is_mat) Y[n0 + seq_len(h), , drop = FALSE] else Y[n0 + seq_len(h)]
      res <- tryCatch({
        fit <- fitfun(Yfit, L = L, r = rr, ...)
        fc <- aeforecast(fit, h = h, method = method, interval = "none")
        c(mse = mean((as.matrix(fc$mean) - as.matrix(Ytest))^2), time = fit$time)
      }, error = function(e) c(mse = NA_real_, time = NA_real_))
      rows[[length(rows) + 1]] <- data.frame(r = rr, origin = n0, mse = res[["mse"]], time = res[["time"]])
      if (verbose) message(sprintf("r = %d, origin %d: MSE %.5g", rr, n0, res[["mse"]]))
    }
  }
  det <- do.call(rbind, rows)
  tab <- stats::aggregate(mse ~ r, data = det, FUN = mean, na.action = stats::na.pass)
  tab$failed <- as.numeric(tapply(is.na(det$mse), det$r, sum)[as.character(tab$r)])
  ok <- is.finite(tab$mse)
  if (!any(ok)) stop("All validation fits failed.")
  out <- list(r = tab$r[ok][which.min(tab$mse[ok])], table = tab, details = det, model = model,
              method = method, h = h, origins = ends, L = L)
  class(out) <- "aessa_select"
  out
}

#' @export
print.aessa_select <- function(x, ...) {
  cat(sprintf("%s latent dimension by rolling-origin validation (h = %d, origins %s, method \"%s\")\n",
              toupper(x$model), x$h, paste(x$origins, collapse = ", "), x$method))
  print(x$table, row.names = FALSE, digits = 5)
  cat(sprintf("selected r = %d\n", x$r))
  invisible(x)
}
