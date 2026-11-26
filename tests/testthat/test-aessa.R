# Tests of the AESSA family: exact reductions of the untrained (linear-start)
# networks to SSA / PrSSA, linear decoder inversion = recurrent SSA forecasting
# and PrSSA window law, training smoke tests, reproducibility, input checks.

skip_if_no_torch <- function() {
  skip_if_not_installed("torch")
  if (!torch::torch_is_installed()) skip("LibTorch is not installed")
}
ae_y <- function(N = 200, sd = 0.7, seed = 1) {
  set.seed(seed)
  tt <- 1:N
  1.5 + 2 * sin(2 * pi * tt / 20) + cos(2 * pi * tt / 7) + rnorm(N, sd = sd)
}

test_that("untrained residual AESSA is Basic SSA of the centred series", {
  skip_if_no_torch()
  skip_if_not_installed("Rssa")
  y <- ae_y()
  L <- 24
  r <- 4
  fit <- aessa(y, L, r, epochs = 0, restarts = 1, valid = 0, seed = 1)
  s <- Rssa::ssa(y - mean(y), L = L, svd.method = "svd")
  rs <- Rssa::reconstruct(s, groups = list(1:r))[[1]]
  rec <- reconstruct(fit)
  expect_lt(max(abs(rec$signal - rs)), 1e-8)
  expect_lt(max(abs(rec$nonlinear)), 1e-14)
  expect_lt(max(abs(rec$linear - rec$signal)), 1e-14)
  expect_equal(as.numeric(rec$mean + rec$signal), as.numeric(fit$signal), tolerance = 1e-12)
  # window loss of the linear start = mean discarded eigenvalue mass per coordinate
  X <- t(Rssa::hankel(y - mean(y), L)) / sqrt(mean((y - mean(y))^2))
  ev <- svd(X)$d^2
  expect_equal(fit$loss, sum(ev[-(1:r)]) / length(X), tolerance = 1e-10)
  expect_equal(fit$ssa_baseline$loss, fit$loss, tolerance = 1e-12)
})

test_that("mean = 'M0' gives uncentred Basic SSA and M2 row-centred SSA", {
  skip_if_no_torch()
  skip_if_not_installed("Rssa")
  y <- ae_y(seed = 2)
  L <- 20
  r <- 3
  f0 <- aessa(y, L, r, mean = "M0", epochs = 0, restarts = 1, valid = 0, seed = 1)
  s0 <- Rssa::ssa(y, L = L, svd.method = "svd")
  expect_lt(max(abs(reconstruct(f0)$signal - Rssa::reconstruct(s0, groups = list(1:r))[[1]])), 1e-8)
  f2 <- aessa(y, L, r, mean = "M2", epochs = 0, restarts = 1, valid = 0, seed = 1)
  # row-centred SSA computed directly: centre the lag coordinates, rank-r SVD, add back, diagonal average
  X <- Rssa::hankel(y, L)
  m <- rowMeans(X)
  sv <- svd(X - m)
  Xr <- sv$u[, 1:r] %*% (sv$d[1:r] * t(sv$v[, 1:r])) + m
  diagavg <- function(M) as.numeric(tapply(M, row(M) + col(M) - 1, mean))
  rec2 <- reconstruct(f2)
  expect_lt(max(abs(as.numeric(rec2$mean + rec2$signal) - as.numeric(f2$signal))), 1e-12)
  expect_lt(max(abs(as.numeric(f2$signal) - diagavg(Xr))), 1e-8)
})

test_that("linear decoder inversion gives minimum-norm recurrent SSA forecasts", {
  skip_if_no_torch()
  skip_if_not_installed("Rssa")
  y <- ae_y(seed = 3)
  L <- 24
  r <- 4
  fit <- aessa(y, L, r, epochs = 0, restarts = 1, valid = 0, seed = 1)
  s <- Rssa::ssa(y - mean(y), L = L, svd.method = "svd")
  fr <- aeforecast(fit, 10, method = "recursive")
  rf <- Rssa::rforecast(s, groups = list(1:r), len = 10, base = "original", only.new = TRUE)
  expect_lt(max(abs(as.numeric(fr$mean) - mean(y) - as.numeric(rf))), 1e-8)
  fw <- aeforecast(fit, 10, method = "window")
  expect_equal(fw$mean[1], fr$mean[1], tolerance = 1e-10)
  expect_true(all(fr$diagnostics$converged))
  expect_equal(fw$method_used, rep("window", 10))
})

test_that("a noise-free finite-rank signal is continued exactly (window and latent)", {
  skip_if_no_torch()
  tt <- 1:150
  sig <- function(t) 2 * sin(2 * pi * t / 12) + 1.2 * cos(2 * pi * t / 30 + 0.4)
  y <- sig(tt)
  fit <- aessa(y, 24, 4, mean = "M0", epochs = 0, restarts = 1, valid = 0, seed = 1)
  truth <- sig(150 + 1:20)
  expect_lt(max(abs(aeforecast(fit, 20, method = "window")$mean - truth)), 1e-7)
  expect_lt(max(abs(aeforecast(fit, 20, method = "latent")$mean - truth)), 1e-7)
})

test_that("untrained VAESSA is PrSSA: sigma2, posterior reconstruction, window law", {
  skip_if_no_torch()
  y <- ae_y(seed = 4)
  L <- 24
  r <- 4
  v <- vaessa(y, L, r, epochs = 0, restarts = 1, valid = 0, seed = 1)
  p <- prssa(y, L, r)
  expect_equal(v$sigma2, p$sigma2, tolerance = 1e-10)
  expect_equal(v$prssa_sigma2, p$sigma2, tolerance = 1e-10)
  expect_lt(max(abs(reconstruct(v)$signal - prreconstruct(p, list(1:r))$F1)), 1e-8)
  pf <- prforecast(p, 12, method = "window", level = c(0.8, 0.95))
  vf <- aeforecast(v, 12, method = "window", interval = "laplace", level = c(0.8, 0.95))
  expect_lt(max(abs(as.numeric(vf$mean) - as.numeric(pf$mean))), 1e-8)
  expect_lt(max(abs(vf$sd - as.numeric(pf$sd))), 1e-8)
  expect_lt(max(abs(as.numeric(vf$lower[["95%"]]) - as.numeric(pf$lower[["95%"]]))), 1e-8)
  expect_lt(max(abs(as.numeric(vf$upper[["80%"]]) - as.numeric(pf$upper[["80%"]]))), 1e-8)
  # the composite ELBO with the exact posterior encoder is the composite likelihood
  # (Monte Carlo over the posterior: agreement up to simulation error)
  expect_lt(abs(v$elbo - p$loglik) / (v$settings$K * v$settings$D), 0.01)
  vm <- aeforecast(v, 6, interval = "mc", npaths = 4000, seed = 1)
  expect_lt(max(abs(vm$sd / vf$sd[1:6] - 1)), 0.08)
})

test_that("untrained block-noise VAESSA is PrMSSA with one noise level per series", {
  skip_if_no_torch()
  y <- ae_y(seed = 5)
  set.seed(6)
  Y <- cbind(a = y, b = 0.5 * y + rnorm(length(y), sd = 1.5))
  v <- vaessa(Y, 20, 3, noise = "block", epochs = 0, restarts = 1, valid = 0, seed = 1)
  p <- prssa(Y, 20, 3, noise = "block")
  expect_equal(unname(v$sigma2), unname(p$sigma2), tolerance = 1e-8)
  pf <- prforecast(p, 6, method = "window")
  vf <- aeforecast(v, 6)
  expect_lt(max(abs(vf$mean - pf$mean)), 1e-8)
  expect_lt(max(abs(vf$sd - pf$sd)), 1e-8)
  expect_equal(dim(vf$mean), c(6L, 2L))
})

test_that("a linear autoencoder trained from a random start recovers SSA", {
  skip_if_no_torch()
  skip_on_cran()
  skip_if_not_installed("Rssa")
  y <- ae_y(seed = 7)
  fit <- aessa(y, 24, 2, arch = aessa_arch(width = 0), init = "random", epochs = 4000, lr = 1e-2,
               patience = 4000, valid = 0, restarts = 1, seed = 3)
  s <- Rssa::ssa(y - mean(y), L = 24, svd.method = "svd")
  rs <- Rssa::reconstruct(s, groups = list(1:2))[[1]]
  expect_lt(max(abs(reconstruct(fit)$signal - rs)), 1e-3)
  expect_lt(fit$loss / fit$ssa_baseline$loss - 1, 1e-6)
})

test_that("AESSA improves on SSA for a stretched square wave at r = 2", {
  skip_if_no_torch()
  skip_on_cran()
  skip_if_not_installed("Rssa")
  set.seed(11)
  N <- 300
  tt <- 1:N
  sig <- tanh(4 * sin(2 * pi * tt / 30))
  x <- sig + rnorm(N, sd = 0.3)
  fa <- aessa(x, 48, 2, seed = 1)
  s <- Rssa::ssa(x - mean(x), L = 48, svd.method = "svd")
  rs <- Rssa::reconstruct(s, groups = list(1:2))[[1]] + mean(x)
  rmse_ssa <- sqrt(mean((rs - sig)^2))
  rmse_ae <- sqrt(mean((fa$signal - sig)^2))
  expect_lt(rmse_ae, 0.7 * rmse_ssa)
  # the validation loss never ends above the SSA start
  expect_true(all(fa$restarts$best_loss <= fa$restarts$start_loss))
  fc <- aeforecast(fa, 24)
  expect_true(all(is.finite(fc$mean)))
  fv <- vaessa(x, 48, 2, seed = 1, restarts = 2)
  expect_lt(sqrt(mean((fv$signal - sig)^2)), 0.7 * rmse_ssa)
  expect_true(all(fv$sigma2 > 0))
  fcv <- aeforecast(fv, 24)
  expect_true(all(fcv$lower[[1]] < fcv$mean & fcv$mean < fcv$upper[[1]]))
})

test_that("fits are reproducible with a seed", {
  skip_if_no_torch()
  y <- ae_y(seed = 8)
  a1 <- aessa(y, 24, 2, seed = 5, restarts = 2, epochs = 60, val_every = 2, patience = 60)
  a2 <- aessa(y, 24, 2, seed = 5, restarts = 2, epochs = 60, val_every = 2, patience = 60)
  expect_equal(a1$par, a2$par, tolerance = 1e-12)
  expect_equal(a1$signal, a2$signal, tolerance = 1e-12)
  v1 <- vaessa(y, 24, 2, seed = 5, restarts = 1, epochs = 40)
  v2 <- vaessa(y, 24, 2, seed = 5, restarts = 1, epochs = 40)
  expect_equal(v1$par, v2$par, tolerance = 1e-12)
})

test_that("matrix and funts inputs keep the input class", {
  skip_if_no_torch()
  y <- ae_y(seed = 9)
  set.seed(9)
  Y <- cbind(u = y, v = -y + rnorm(length(y)))
  fm <- aessa(Y, 20, 2, seed = 1, restarts = 1, epochs = 30)
  expect_equal(dim(fm$signal), dim(Y))
  expect_equal(dim(aeforecast(fm, 5)$mean), c(5L, 2L))
  grid <- seq(0, 1, length.out = 25)
  M <- sapply(1:50, function(i) sin(2 * pi * grid) * cos(2 * pi * i / 12) + rnorm(25, sd = 0.2))
  Yf <- funts(X = M, basisobj = fda::create.bspline.basis(c(0, 1), 9))
  ff <- vaessa(Yf, 12, 2, seed = 1, restarts = 1, epochs = 0, valid = 0)
  expect_s3_class(ff$signal, "funts")
  expect_equal(ff$sigma2, prssa(Yf, 12, 2)$sigma2, tolerance = 1e-8)
  expect_s3_class(aeforecast(ff, 2)$mean, "funts")
})

test_that("methods run and inputs are validated", {
  skip_if_no_torch()
  y <- ae_y(seed = 10)
  fit <- aessa(y, 24, 2, seed = 1, restarts = 1, epochs = 20)
  expect_output(print(fit), "AESSA fit")
  expect_output(print(summary(fit)), "restarts")
  pdf(NULL)
  plot(fit, "loss")
  plot(fit, "latent")
  plot(fit, "components")
  grDevices::dev.off()
  expect_s3_class(predict(fit, h = 3), "aeforecast")
  expect_error(aessa(y, L = 1, r = 2), "L")
  expect_error(aessa(y, 24, 0), "r")
  expect_error(aessa(y, 24, 2, valid = 0.7), "valid")
  expect_error(aessa(c(y, NA), 24, 2), "missing")
  expect_error(aessa_arch(width = 0, residual = FALSE), "residual")
  expect_error(aessa(y, 24, 2, arch = aessa_arch(residual = FALSE)), "residual")
  expect_error(aeforecast(fit, 5, interval = "laplace"), "VAESSA")
  expect_error(aeforecast(fit, 0), "h")
  expect_error(aessa(y, 24, 2, arch = list(width = 3)), "aessa_arch")
})
