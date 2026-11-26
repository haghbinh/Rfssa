# Tests of the PrSSA family: exact reductions to SSA / MSSA / FSSA, ECME
# monotonicity, forecasting limits, the sandwich penalty and basic coverage.

sim_y <- function(N = 200, sd = 0.7, seed = 1) {
  set.seed(seed)
  tt <- 1:N
  1.5 + 2 * sin(2 * pi * tt / 20) + cos(2 * pi * tt / 7) + rnorm(N, sd = sd)
}
proj_dist <- function(A, B) {
  qa <- qr.Q(qr(A))
  qb <- qr.Q(qr(B))
  norm(tcrossprod(qa) - tcrossprod(qb), "F")
}
hankel_traj <- function(y, L) outer(seq_len(L), seq_len(length(y) - L + 1), function(i, j) y[i + j - 1])

test_that("isotropic PrSSA with mean = 'zero' is Basic SSA (subspace, sigma2)", {
  skip_if_not_installed("Rssa")
  y <- sim_y()
  N <- length(y); L <- 24; r <- 4; K <- N - L + 1
  fit <- prssa(y, L, r, mean = "zero")
  s <- Rssa::ssa(y, L = L, svd.method = "svd")
  expect_lt(proj_dist(fit$W, s$U[, 1:r]), 1e-8)
  w <- pmin(seq_len(N), L, K, N:1)
  sig2 <- (sum(w * y^2) - sum(s$sigma[1:r]^2)) / (K * (L - r))
  expect_equal(fit$sigma2, sig2, tolerance = 1e-10)
  ev <- eigen(tcrossprod(hankel_traj(y, L)) / K, symmetric = TRUE)$values
  expect_equal(fit$sigma2, mean(ev[-(1:r)]), tolerance = 1e-10)
  expect_equal(fit$eigenvalues[1:r] * K, s$sigma[1:r]^2, tolerance = 1e-10)
})

test_that("projection reconstruction is SSA grouping; posterior is shrunken SSA", {
  skip_if_not_installed("Rssa")
  y <- sim_y()
  L <- 24; r <- 4
  fit <- prssa(y, L, r)
  s <- Rssa::ssa(y - mean(y), L = L, svd.method = "svd")
  rp <- prreconstruct(fit, list(1:2, 3:4), type = "projection")
  rs <- Rssa::reconstruct(s, groups = list(1:2, 3:4))
  expect_lt(max(abs(rp$F1 - rs$F1)), 1e-8)
  expect_lt(max(abs(rp$F2 - rs$F2)), 1e-8)
  expect_equal(as.numeric(rp$mean), rep(mean(y), length(y)))
  el <- Rssa::reconstruct(s, groups = as.list(1:r))
  shr <- Reduce(`+`, Map(function(w, comp) w * comp, fit$omega, el[1:r]))
  expect_lt(max(abs(prreconstruct(fit)$F1 - shr)), 1e-8)
  expect_true(all(fit$omega > 0 & fit$omega < 1))
})

test_that("block PrMSSA: ECME is monotone and its fixed point is weighted MSSA", {
  skip_if_not_installed("Rssa")
  set.seed(2)
  N <- 300; L <- 40; K <- N - L + 1; tt <- 1:N
  y1 <- 2 * sin(2 * pi * tt / 20) + cos(2 * pi * tt / 7) + rnorm(N, sd = 0.3)
  y2 <- 1.5 * sin(2 * pi * tt / 20 + 0.5) + 0.8 * cos(2 * pi * tt / 7 + 1) + rnorm(N, sd = 2)
  fit <- prssa(cbind(y1, y2), L, 4, noise = "block", control = list(tol = 1e-12, trace = TRUE))
  expect_true(fit$converged)
  expect_gte(min(diff(fit$trace)), -1e-8)
  expect_equal(unname(sqrt(fit$sigma2)), c(0.3, 2), tolerance = 0.15)
  # fixed point sigma_b^2 = RSS_b / (K (D_b - h_b))
  X <- .pr_fit_trajectory(fit)
  Xt <- X / sqrt(fit$psi)
  blk <- rep(1:2, each = L)
  RSS <- tapply(rowSums((X - X %*% tcrossprod(fit$V))^2), blk, sum)
  hb <- tapply(rowSums(fit$U^2), blk, sum)
  expect_equal(as.numeric(fit$sigma2), as.numeric(RSS / (K * (L - hb))), tolerance = 1e-8)
  # vertical PrMSSA(L) = Rssa horizontal MSSA with window K on the scaled series
  sc <- sqrt(fit$sigma2)
  mR <- Rssa::ssa(list((y1 - mean(y1)) / sc[1], (y2 - mean(y2)) / sc[2]), L = K, kind = "mssa", svd.method = "svd")
  expect_equal(fit$s, mR$sigma[1:4], tolerance = 1e-8)
  rR <- Rssa::reconstruct(mR, groups = list(1:4))[[1]]
  rp <- prreconstruct(fit, list(1:4), type = "projection")$F1
  expect_lt(max(abs(rp[, 1] - as.numeric(rR[[1]]) * sc[1])), 1e-8)
  expect_lt(max(abs(rp[, 2] - as.numeric(rR[[2]]) * sc[2])), 1e-8)
})

test_that("sigma -> 0 forecast limits are the SSA LRF, rforecast(base = 'original') and vforecast", {
  skip_if_not_installed("Rssa")
  y <- sim_y()
  L <- 24; r <- 4
  fit <- prssa(y, L, r, mean = "zero")
  s <- Rssa::ssa(y, L = L, svd.method = "svd")
  lim <- .pr_noise_limit(fit, 1e-12)
  os <- .pr_fc_onestep(lim$pred$W, lim$pred$psi, .pr_record_rows(1L, L), L)
  expect_lt(max(abs(as.numeric(os$WL %*% os$B) - as.numeric(Rssa::lrr(s, groups = list(1:r))))), 1e-6)
  f1 <- prforecast(lim, 1, method = "recursive", bands = "pointwise")
  rf <- Rssa::rforecast(s, groups = list(1:r), len = 1, base = "original", only.new = TRUE)
  expect_lt(abs(as.numeric(f1$mean) - as.numeric(rf)), 1e-6)
  fv <- prforecast(lim, 8, method = "vector")
  vf <- Rssa::vforecast(s, groups = list(1:r), len = 8, only.new = TRUE)
  expect_lt(max(abs(as.numeric(fv$mean) - as.numeric(vf))), 1e-6)
  # at the fitted noise level the one-step mean is a ridge LRF (closed form)
  U <- fit$U; om <- fit$omega
  a_pr <- U[-L, ] %*% (om * U[L, ]) / (1 - sum(om * U[L, ]^2))
  fr <- prforecast(fit, 1, method = "recursive", bands = "pointwise")
  expect_equal(as.numeric(fr$mean), sum(a_pr * tail(y, L - 1)), tolerance = 1e-10)
  expect_equal(fr$sd^2, fit$sigma2 / (1 - sum(om * U[L, ]^2)), tolerance = 1e-10)
})

test_that("vector mode shrinks once: exact V-forecast of span(W) and no amplitude damping", {
  # Golyandina's V-forecast operator from a basis U, seeded with trajectory columns Z
  gvf <- function(U, Z, L, h) {
    Un <- U[-L, , drop = FALSE]; p <- U[L, ]; nu2 <- sum(p^2)
    R <- Un %*% p / (1 - nu2)
    Pi <- tcrossprod(Un) + (1 - nu2) * tcrossprod(R)
    K <- ncol(Z); last <- Z[, K]
    for (i in seq_len(h + L - 1)) { yb <- last[2:L]; last <- c(Pi %*% yb, sum(R * yb)); Z <- cbind(Z, last) }
    sapply(K + L - 1 + seq_len(h), function(q) {
      j <- max(1, q - L + 1):min(ncol(Z), q); mean(Z[cbind(q - j + 1, j)]) })
  }
  sig <- function(t) 2 * sin(2 * pi * t / 20) + cos(2 * pi * t / 7)
  set.seed(2)
  y <- sig(1:300) + rnorm(300, sd = 1.5)
  fit <- prssa(y, 40, 4, mean = "zero")
  Shat <- (sqrt(fit$psi) * fit$U) %*% ((fit$omega * fit$s) * t(fit$V))
  expect_lt(max(abs(as.numeric(prforecast(fit, 15, method = "vector")$mean) - gvf(fit$U, Shat, 40, 15))), 1e-8)
  # long harmonic series: forecast amplitude is not systematically damped
  N <- 1000; L <- 120; h <- 120
  ratio <- sapply(1:3, function(seed) {
    set.seed(seed)
    y <- sig(1:N) + rnorm(N, sd = 1.5)
    f <- as.numeric(prforecast(prssa(y, L, 4), h, method = "vector")$mean)
    sqrt(mean((f - mean(y))^2) / mean(sig(N + 1:h)^2))
  })
  expect_gt(mean(ratio), 0.8)
  expect_lt(mean(ratio), 1.2)
})

test_that("marginal window laws: h = 1 unchanged, each horizon equals a separate joint call", {
  y <- sim_y()
  L <- 24; r <- 4; N <- length(y)
  fit <- prssa(y, L, r)
  H <- 6
  fm <- prforecast(fit, H, bands = "pointwise")
  expect_identical(prforecast(fit, 1, bands = "pointwise")$mean,
                   prforecast(fit, 1, joint = TRUE, bands = "pointwise")$mean)
  expect_identical(prforecast(fit, 1, bands = "pointwise")$sd,
                   prforecast(fit, 1, joint = TRUE, bands = "pointwise")$sd)
  for (k in seq_len(H)) {
    fk <- prforecast(fit, k, joint = TRUE, bands = "pointwise")
    expect_equal(as.numeric(fm$mean)[k], as.numeric(fk$mean)[k], tolerance = 1e-12)
    expect_equal(fm$sd[k], fk$sd[k], tolerance = 1e-12)
  }
  # horizon k: Gaussian conditioning of the window ending at N + k on its first L - k lags
  C <- tcrossprod(fit$W) + diag(fit$sigma2, L)
  yc <- y - mean(y)
  for (k in c(1, 3, 6)) {
    P <- 1:(L - k)
    B <- C[L, P, drop = FALSE] %*% solve(C[P, P])
    expect_equal(as.numeric(fm$mean)[k], as.numeric(B %*% yc[(N - L + k + 1):N]) + mean(y), tolerance = 1e-10)
    expect_equal(fm$sd[k]^2, as.numeric(C[L, L] - B %*% C[P, L]), tolerance = 1e-10)
  }
  # several series: per-horizon equality and the joint cross-horizon covariance
  set.seed(5)
  fb <- prssa(cbind(a = y, b = 0.5 * y + rnorm(N, sd = 0.8)), L, r, noise = "block")
  gm <- prforecast(fb, 4, bands = "pointwise")
  g3 <- prforecast(fb, 3, joint = TRUE, bands = "pointwise")
  gj <- prforecast(fb, 4, joint = TRUE, bands = "pointwise")
  expect_equal(gm$mean[3, ], g3$mean[3, ], tolerance = 1e-12)
  expect_equal(gm$cov[[3]], g3$cov[[3]], tolerance = 1e-12)
  expect_equal(gj$joint_cov[7:8, 7:8], gj$cov[[4]], tolerance = 1e-12)
  expect_null(gm$joint_cov)
})

test_that("window-conditional and recursive laws agree with Gaussian conditioning", {
  y <- sim_y()
  L <- 24; r <- 4; N <- length(y)
  fit <- prssa(y, L, r)
  C <- tcrossprod(fit$W) + diag(fit$sigma2, L)
  h <- 3
  P <- 1:(L - h); Fi <- (L - h + 1):L
  xP <- (y - mean(y))[(N - L + h + 1):N]
  B <- C[Fi, P] %*% solve(C[P, P])
  fw <- prforecast(fit, h, method = "window", bands = "pointwise", joint = TRUE)
  expect_equal(as.numeric(fw$mean), as.numeric(B %*% xP) + mean(y), tolerance = 1e-10)
  expect_equal(fw$sd^2, diag(C[Fi, Fi] - B %*% C[P, Fi]), tolerance = 1e-10)
  # recursive variances = AR(L-1) impulse-response formula
  a <- as.numeric(prforecast(fit, 1, "recursive", bands = "pointwise")$sd)^2
  fr <- prforecast(fit, 6, method = "recursive", bands = "pointwise")
  U <- fit$U; om <- fit$omega
  ar <- as.numeric(U[-L, ] %*% (om * U[L, ]) / (1 - sum(om * U[L, ]^2)))
  psiw <- numeric(6); psiw[1] <- 1
  for (i in 2:6) { k <- i - 1; psiw[i] <- sum(rev(ar)[1:min(k, L - 1)] * psiw[k:max(1, k - L + 2)]) }
  expect_equal(fr$sd^2, a * cumsum(psiw^2), tolerance = 1e-10)
  expect_equal(prforecast(fit, 30)$method_used, rep(c("window", "recursive"), c(L - r, 30 - (L - r))))
})

test_that("reduced sandwich penalty equals the explicit tr(J H^+)", {
  set.seed(3)
  N <- 120; L <- 12; K <- N - L + 1; tt <- 1:N
  y1 <- sin(2 * pi * tt / 12) + rnorm(N, sd = 0.3)
  y2 <- 0.8 * sin(2 * pi * tt / 12 + 1) + rnorm(N, sd = 0.8)
  fit <- prssa(cbind(y1, y2), L, 2, noise = "block", control = list(tol = 1e-12))
  X <- .pr_fit_trajectory(fit); W <- fit$W; psi <- fit$psi; blk <- rep(1:2, each = L)
  D <- nrow(X); r <- ncol(W)
  Ci <- solve(tcrossprod(W) + diag(psi)); P <- Ci %*% W; CiE <- Ci %*% X; Zw <- crossprod(X, P)
  SW <- do.call(cbind, lapply(1:r, function(a) t(CiE) * Zw[, a])) - rep(as.numeric(P), each = K)
  Ss <- sapply(1:2, function(b) 0.5 * (colSums(CiE[blk == b, , drop = FALSE]^2) - sum(diag(Ci)[blk == b])))
  Sm <- cbind(SW, Ss)
  Q <- crossprod(W, P); iv <- rep(1:D, r); bv <- rep(1:r, each = D); F1 <- P[iv, bv]
  HWW <- F1 * t(F1) + kronecker(Q, Ci)
  HWs <- sapply(1:2, function(b) as.numeric(Ci[, blk == b] %*% P[blk == b, ]))
  Hss <- outer(1:2, 1:2, Vectorize(function(b, c) 0.5 * sum(Ci[blk == b, blk == c]^2)))
  Hp <- .pr_pinv(rbind(cbind(HWW, HWs), cbind(t(HWs), Hss)))
  pe_explicit <- .pr_hac_trace(Sm, Sm %*% Hp, .pr_hac_weights(L, K))
  expect_equal(.pr_peff(fit), pe_explicit, tolerance = 1e-8)
})

test_that("prrank selects the true rank on a clear example", {
  y <- sim_y(N = 300, sd = 0.5)
  rk <- prrank(y, 24, rmax = 6)
  expect_equal(rk$selected, 4)
  expect_equal(prrank(y, 24, rmax = 6, criterion = "clbic_nominal")$selected, 4)
})

test_that("one-step prediction intervals have roughly nominal coverage", {
  set.seed(4)
  hit <- logical(60)
  for (i in seq_along(hit)) {
    tt <- 1:151
    sig <- 2 * sin(2 * pi * tt / 20 + runif(1, 0, 2 * pi)) + cos(2 * pi * tt / 7 + runif(1, 0, 2 * pi))
    y <- sig + rnorm(151, sd = 0.6)
    fc <- prforecast(prssa(y[1:150], 24, 4), 1, level = 0.9, bands = "pointwise")
    hit[i] <- y[151] >= fc$lower[[1]] && y[151] <= fc$upper[[1]]
  }
  expect_gte(mean(hit), 0.75)
})

test_that("PrFSSA on an orthonormal Fourier basis reproduces Rfssa::fssa eigenvalues", {
  set.seed(5)
  u <- seq(0, 1, length.out = 41); N <- 100; L <- 10; K <- N - L + 1
  Xf <- outer(sin(2 * pi * u), sin(2 * pi * (1:N) / 10)) + outer(cos(2 * pi * u), cos(2 * pi * (1:N) / 10)) +
    matrix(rnorm(41 * N, sd = 0.3), 41)
  Y <- funts(Xf, fda::create.fourier.basis(c(0, 1), 5))
  fit <- prssa(Y, L, 3, mean = "zero")
  utils::capture.output(fs <- fssa(Y, L = L, ntriples = 3))
  expect_equal(fit$eigenvalues[1:3] * K, fs$values[1:3], tolerance = 1e-6)
})

test_that("functional input returns funts objects; hetero and smoothing paths work", {
  set.seed(6)
  u <- seq(0, 1, length.out = 30); N <- 120; tt <- 1:N
  sdv <- 0.1 + 0.5 * u
  Xf <- outer(exp(u^2), cos(2 * pi * tt / 10)) + outer(cos(4 * pi * u), sin(2 * pi * tt / 10)) +
    sdv * matrix(rnorm(30 * N), 30)
  Y <- funts(Xf, fda::create.bspline.basis(c(0, 1), 12))
  fit <- prssa(Y, L = 10, r = 2)
  rec <- prreconstruct(fit)
  expect_s3_class(rec$F1, "funts")
  expect_equal(rec$F1$N, N)
  fc <- prforecast(fit, h = 3, nsim = 200)
  expect_s3_class(fc$mean, "funts")
  expect_s3_class(fc$intervals$simultaneous$upper[[2]], "funts")
  expect_true(all(fc$quantile[2, ] >= stats::qnorm(0.975)))
  fh <- prssa(Y, L = 10, r = 2, noise = "hetero", hetero_df = 6, control = list(trace = TRUE))
  expect_gte(min(diff(fh$trace)), -1e-8)
  s_hat <- sqrt(fh$sigma2[[1]])
  expect_gt(stats::cor(s_hat, sdv), 0.9)
  fs <- prssa(Y, L = 10, r = 2, alpha = 1e-13)
  expect_lt(proj_dist(fs$W, fit$W), 1e-6)
  fs2 <- prssa(Y, L = 10, r = 2, alpha = 1e-4)
  expect_true(fs2$stage2$sigma2 >= fs2$sigma2)
  expect_s3_class(prforecast(fs2, 2, bands = "pointwise")$mean, "funts")
  # hybrid list: one funts variable and one scalar series
  z <- sin(2 * pi * tt / 10) + rnorm(N, sd = 0.2)
  fhy <- prssa(list(curves = Y, z = z), L = 10, r = 2, noise = "block")
  expect_length(fhy$sigma2, 2)
  rh <- prreconstruct(fhy)$F1
  expect_s3_class(rh$curves, "funts")
  expect_length(rh$z, N)
})
