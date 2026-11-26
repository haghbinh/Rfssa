# =============================================================================
# Internal utilities for the PrSSA family (probabilistic singular spectrum
# analysis): data preparation, lag embedding, diagonal averaging, the noise
# structure and conversion of results back to the class of the input.
#
# Records are stacked vertically (VMSSA / MFSSA layout): block-outer, then the
# basis (or grid) coordinate q of the block, then the lag l (inner).  The
# coordinate (b, q, l) of a record therefore has row index
#   offset_b + (q - 1) * L + l,   offset_b = L * sum_{c < b} d_c.
# Functional blocks are whitened with the symmetric root of the exact L2 Gram
# matrix of their basis, so Euclidean geometry of the working coordinates is
# the L2 geometry of the curves.
# =============================================================================

.pr_control <- function(control) {
  ctrl <- list(maxit = 500, tol = 1e-8, trace = FALSE, floor = 1e-6,
               nextra = 10, svd = "auto", svds_min = 300)
  if (!is.null(control)) {
    if (!is.list(control)) stop("`control` must be a list.")
    ctrl[names(control)] <- control
  }
  ctrl
}

# lag-embedding index: L x K matrix whose column j holds j, ..., j + L - 1
.pr_hankel_index <- function(N, L) {
  outer(seq_len(L), seq_len(N - L + 1L) - 1L, "+")
}

# number of windows that contain time t = 1, ..., N
.pr_hankel_weights <- function(N, L) {
  K <- N - L + 1L
  pmin(seq_len(N), L, K, N:1)
}

# stacked (vertical) trajectory matrix from a list of d_b x N working series
.pr_trajectory <- function(Zlist, L) {
  N <- ncol(Zlist[[1]])
  idx <- .pr_hankel_index(N, L)
  do.call(rbind, lapply(Zlist, function(Z) {
    do.call(rbind, lapply(seq_len(nrow(Z)), function(q) {
      z <- Z[q, ]
      matrix(z[idx], nrow = L)
    }))
  }))
}

# diagonal averaging of a stacked D x K' trajectory: list of d_b x (K'+L-1)
.pr_hankelize <- function(Tm, L, dims) {
  Kc <- ncol(Tm)
  N <- Kc + L - 1L
  w <- .pr_hankel_weights(N, L)
  out <- vector("list", length(dims))
  off <- 0L
  for (b in seq_along(dims)) {
    d <- dims[b]
    S <- matrix(0, d, N)
    for (l in seq_len(L)) {
      rows <- off + (seq_len(d) - 1L) * L + l
      cols <- l:(l + Kc - 1L)
      S[, cols] <- S[, cols, drop = FALSE] + Tm[rows, , drop = FALSE]
    }
    out[[b]] <- sweep(S, 2, w, "/")
    off <- off + d * L
  }
  out
}

# apply the block-diagonal map diag_b(A_b (x) I_L) to the rows of a D x n matrix
.pr_apply_blockdiag <- function(Xm, Alist, L, dims) {
  out <- Xm
  off <- 0L
  for (b in seq_along(dims)) {
    d <- dims[b]
    A <- Alist[[b]]
    if (!is.null(A) && !(d == 1L && A[1, 1] == 1)) {
      for (l in seq_len(L)) {
        rows <- off + (seq_len(d) - 1L) * L + l
        out[rows, ] <- A %*% Xm[rows, , drop = FALSE]
      }
    }
    off <- off + d * L
  }
  out
}

# record row index of (time-vector position p, lag l): d_tot x L matrix
.pr_record_rows <- function(dims, L) {
  dt <- sum(dims)
  ri <- matrix(0L, dt, L)
  p <- 0L
  off <- 0L
  for (b in seq_along(dims)) {
    for (q in seq_len(dims[b])) {
      p <- p + 1L
      ri[p, ] <- off + (q - 1L) * L + seq_len(L)
    }
    off <- off + dims[b] * L
  }
  ri
}

# symmetric square root and inverse square root of a symmetric PSD matrix
.pr_msqrt <- function(A, inverse = FALSE) {
  e <- eigen((A + t(A)) / 2, symmetric = TRUE)
  v <- pmax(e$values, 0)
  if (inverse) {
    if (any(v <= .Machine$double.eps * max(v) * 10)) {
      stop("The Gram matrix of the basis is (numerically) singular.")
    }
    v <- 1 / sqrt(v)
  } else {
    v <- sqrt(v)
  }
  e$vectors %*% (v * t(e$vectors))
}

# truncated SVD: base svd for small problems, RSpectra::svds for large ones
.pr_svd <- function(X, k, ctrl) {
  m <- min(dim(X))
  k <- min(k, m)
  use_svds <- identical(ctrl$svd, "svds") ||
    (identical(ctrl$svd, "auto") && m > ctrl$svds_min && k < m / 3)
  if (use_svds) {
    s <- RSpectra::svds(X, k = k, nu = k, nv = k)
    return(list(d = s$d, u = s$u, v = s$v))
  }
  s <- svd(X, nu = k, nv = k)
  list(d = s$d[seq_len(k)], u = s$u, v = s$v)
}

# orthonormal basis of the column space of A
.pr_orth <- function(A, tol = 1e-10) {
  if (ncol(A) == 0) return(A)
  s <- svd(A)
  keep <- s$d > tol * max(s$d)
  s$u[, keep, drop = FALSE]
}

# ------------------------------------------------------------------------------
# Input preparation
# ------------------------------------------------------------------------------

# exact L2 Gram matrix of an fda basis (1D) or of a tensor basis (2D)
.pr_gram <- function(basis, dimS) {
  if (dimS == 1) {
    if (!fda::is.basis(basis)) {
      stop("PrFSSA needs `funts` variables built on an fda basis (empirical bases are not supported).")
    }
    G <- fda::inprod(basis, basis)
  } else {
    if (!(is.list(basis) && length(basis) == 2 && all(vapply(basis, fda::is.basis, logical(1))))) {
      stop("PrFSSA needs two-dimensional `funts` variables built on two fda bases.")
    }
    G <- kronecker(fda::inprod(basis[[1]], basis[[1]]), fda::inprod(basis[[2]], basis[[2]]))
  }
  G <- as.matrix(G)
  (G + t(G)) / 2
}

.pr_scalar_block <- function(y, name) {
  list(kind = "scalar", grid = FALSE, name = name, d = 1L,
       Z = matrix(as.numeric(y), nrow = 1L),
       R = matrix(1), Rinv = matrix(1), S = NULL)
}

.pr_funts_blocks <- function(Y, noise, alpha, element = 1L) {
  p <- length(Y$dimSupp)
  vn <- Y$vnames
  out <- vector("list", p)
  for (j in seq_len(p)) {
    nm <- if (!is.null(vn) && length(vn) >= j && !is.na(vn[j])) vn[j] else paste("Variable", j)
    dimS <- Y$dimSupp[[j]]
    C <- as.matrix(Y$coefs[[j]])
    if (identical(noise, "hetero")) {
      if (dimS != 1) stop("noise = \"hetero\" needs one-dimensional functional variables.")
      if (alpha > 0) stop("Smoothing (alpha > 0) is not available with noise = \"hetero\".")
      grid <- as.numeric(Y$argval[[j]])
      Zg <- Y$B_mat[[j]] %*% C
      out[[j]] <- list(kind = "functional", grid = TRUE, name = nm, d = length(grid),
                       Z = Zg, R = NULL, Rinv = NULL, S = NULL, s = grid,
                       var = j, element = element)
    } else {
      G <- .pr_gram(Y$basis[[j]], dimS)
      Rt <- .pr_msqrt(G)
      Rinv <- .pr_msqrt(G, inverse = TRUE)
      S <- NULL
      if (alpha > 0) {
        if (dimS != 1) stop("Smoothing (alpha > 0) is available for one-dimensional functional variables only.")
        Dn <- as.matrix(fda::inprod(Y$basis[[j]], Y$basis[[j]], 2, 2))
        Om <- Rinv %*% Dn %*% Rinv
        e <- eigen((Om + t(Om)) / 2, symmetric = TRUE)
        S <- e$vectors %*% ((1 / sqrt(1 + alpha * pmax(e$values, 0))) * t(e$vectors))
      }
      out[[j]] <- list(kind = "functional", grid = FALSE, name = nm, d = nrow(C),
                       Z = Rt %*% C, R = Rt, Rinv = Rinv, S = S,
                       var = j, element = element)
    }
  }
  out
}

.pr_prepare <- function(Y, noise, alpha, grid = NULL) {
  meta <- list()
  if (inherits(Y, "funts")) {
    type <- "funts"
    blocks <- .pr_funts_blocks(Y, noise, alpha)
    meta$elements <- list(list(kind = "funts", obj = Y, blocks = seq_along(blocks)))
  } else if (is.list(Y) && !is.data.frame(Y)) {
    type <- "list"
    if (length(Y) == 0) stop("`Y` is an empty list.")
    nms <- names(Y)
    blocks <- list()
    meta$elements <- vector("list", length(Y))
    for (i in seq_along(Y)) {
      nm <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else paste("Series", i)
      yi <- Y[[i]]
      if (inherits(yi, "funts")) {
        bl <- .pr_funts_blocks(yi, noise, alpha, element = i)
        meta$elements[[i]] <- list(kind = "funts", obj = yi, name = nm,
                                   blocks = length(blocks) + seq_along(bl))
        blocks <- c(blocks, bl)
      } else if (is.numeric(yi) && is.null(dim(yi))) {
        blocks <- c(blocks, list(.pr_scalar_block(yi, nm)))
        meta$elements[[i]] <- list(kind = "scalar", name = nm, blocks = length(blocks))
      } else {
        stop("Elements of a list `Y` must be numeric vectors or `funts` objects.")
      }
    }
  } else if (is.matrix(Y) || is.data.frame(Y)) {
    type <- "matrix"
    Ym <- as.matrix(Y)
    if (!is.numeric(Ym)) stop("`Y` must be numeric.")
    cn <- colnames(Ym)
    if (is.null(cn)) cn <- paste("Series", seq_len(ncol(Ym)))
    blocks <- lapply(seq_len(ncol(Ym)), function(j) .pr_scalar_block(Ym[, j], cn[j]))
    meta$colnames <- cn
    meta$tsp <- if (stats::is.ts(Y)) stats::tsp(Y) else NULL
    if (identical(noise, "hetero")) {
      # one curve per row: the columns are the grid points of a single functional variable
      if (alpha > 0) stop("Smoothing (alpha > 0) is not available with noise = \"hetero\".")
      m <- ncol(Ym)
      if (is.null(grid)) grid <- seq(0, 1, length.out = m)
      grid <- as.numeric(grid)
      if (length(grid) != m || any(!is.finite(grid)) || (m > 1 && any(diff(grid) <= 0))) {
        stop("`grid` must be an increasing numeric vector with one value per column of `Y`.")
      }
      meta$grid <- grid
    }
  } else if (is.numeric(Y)) {
    type <- "vector"
    blocks <- list(.pr_scalar_block(Y, "y"))
    meta$tsp <- if (stats::is.ts(Y)) stats::tsp(Y) else NULL
  } else {
    stop("`Y` must be a numeric vector or ts, a numeric matrix, a list of series, or a `funts` object.")
  }
  Ns <- vapply(blocks, function(b) ncol(b$Z), integer(1))
  if (length(unique(Ns)) != 1) stop("All series in `Y` must have the same length.")
  if (any(!is.finite(unlist(lapply(blocks, function(b) b$Z))))) {
    stop("`Y` contains missing or non-finite values.")
  }
  if (identical(noise, "hetero") && type != "matrix" &&
      !any(vapply(blocks, function(b) isTRUE(b$grid), logical(1)))) {
    stop("noise = \"hetero\" needs a functional (`funts`) variable or a matrix of curves on a grid.")
  }
  list(type = type, blocks = blocks, N = Ns[1], meta = meta,
       dims = vapply(blocks, function(b) as.integer(b$d), integer(1)),
       names = vapply(blocks, function(b) b$name, character(1)))
}

# ------------------------------------------------------------------------------
# Noise structure: every coordinate belongs to a noise unit; units are grouped
# into "const" groups (one variance) or "loglin" groups (log variance linear in
# a B-spline basis of the grid argument, heteroscedastic PrFSSA).
# ------------------------------------------------------------------------------
.pr_hetero_basis <- function(s, q) {
  if (q <= 1) return(matrix(1, length(s), 1))
  if (q < 4) return(cbind(1, stats::poly(s, degree = q - 1)))
  a <- min(s)
  b <- max(s)
  inner <- if (q > 4) seq(a, b, length.out = q - 2)[-c(1, q - 2)] else numeric(0)
  splines::splineDesign(knots = c(rep(a, 4), inner, rep(b, 4)), x = s, ord = 4, outer.ok = TRUE)
}

.pr_noise_structure <- function(dat, noise, L, K, hetero_df, hetero_lambda) {
  dims <- dat$dims
  D <- sum(dims) * L
  coord_block <- rep(seq_along(dims), dims * L)
  groups <- list()
  if (noise == "iso") {
    unit <- rep(1L, D)
    groups[[1]] <- list(type = "const", units = 1L, name = "all")
  } else if (noise == "block") {
    unit <- coord_block
    for (b in seq_along(dims)) groups[[b]] <- list(type = "const", units = b, name = dat$names[b])
  } else if (dat$type == "matrix") {
    # curves on a grid: one noise unit per column (grid point), tied over the lags,
    # log sigma^2(s_j) = b(s_j)' eta
    m <- length(dims)
    unit <- coord_block
    s <- dat$meta$grid
    q <- max(1L, min(as.integer(hetero_df), m))
    Bd <- .pr_hetero_basis(s, q)
    Pen <- if (ncol(Bd) >= 3) crossprod(diff(diag(ncol(Bd)), differences = 2)) else matrix(0, ncol(Bd), ncol(Bd))
    lam <- if (is.null(hetero_lambda)) 1 else hetero_lambda
    groups[[1]] <- list(type = "loglin", units = seq_len(m), Bd = Bd, P = Pen,
                        lambda = lam * K * L / 2, name = "Y", block = NA_integer_, s = s)
  } else {
    unit <- integer(D)
    u <- 0L
    off <- 0L
    for (b in seq_along(dims)) {
      blk <- dat$blocks[[b]]
      if (isTRUE(blk$grid)) {
        ids <- u + seq_len(dims[b])
        unit[off + seq_len(dims[b] * L)] <- rep(ids, each = L)
        q <- max(1L, min(as.integer(hetero_df), dims[b]))
        Bd <- .pr_hetero_basis(blk$s, q)
        Pen <- if (ncol(Bd) >= 3) crossprod(diff(diag(ncol(Bd)), differences = 2)) else matrix(0, ncol(Bd), ncol(Bd))
        lam <- if (is.null(hetero_lambda)) 1 else hetero_lambda
        groups[[length(groups) + 1]] <- list(type = "loglin", units = ids, Bd = Bd, P = Pen,
                                             lambda = lam * K * L / 2, name = blk$name, block = b)
        u <- u + dims[b]
      } else {
        u <- u + 1L
        unit[off + seq_len(dims[b] * L)] <- u
        groups[[length(groups) + 1]] <- list(type = "const", units = u, name = blk$name, block = b)
      }
      off <- off + dims[b] * L
    }
  }
  unit_n <- tabulate(unit, max(unit))
  list(unit = unit, unit_n = unit_n, groups = groups, D = D,
       iso = (length(groups) == 1 && groups[[1]]$type == "const"))
}

# number of free noise parameters
.pr_noise_npar <- function(nz) {
  sum(vapply(nz$groups, function(g) if (g$type == "const") 1 else ncol(g$Bd), numeric(1)))
}

# ------------------------------------------------------------------------------
# Output conversion
# ------------------------------------------------------------------------------

# time stamps of a funts object, extended beyond N when needed
.pr_funts_times <- function(obj, idx) {
  tm <- obj$time
  N <- length(tm)
  if (all(idx <= N)) return(tm[idx])
  dtm <- if (N > 1) tm[N] - tm[N - 1] else 1
  out <- tm[pmin(idx, N)]
  beyond <- idx > N
  out[beyond] <- tm[N] + (idx[beyond] - N) * dtm
  out
}

# build a funts object from per-variable series in working coordinates
.pr_make_funts <- function(obj, blocks, series, idx) {
  p <- length(blocks)
  grid_mode <- isTRUE(blocks[[1]]$grid)
  Xl <- vector("list", p)
  argv <- vector("list", p)
  for (j in seq_len(p)) {
    blk <- blocks[[j]]
    v <- blk$var
    if (grid_mode) {
      Xl[[j]] <- series[[j]]
      argv[[j]] <- blk$s
    } else {
      Xl[[j]] <- blk$Rinv %*% series[[j]]
      if (obj$dimSupp[[v]] == 1) {
        argv[[j]] <- as.numeric(obj$argval[[v]])
      } else {
        argv[[j]] <- list(unique(obj$argval[[v]][, 1]), unique(obj$argval[[v]][, 2]))
      }
    }
  }
  vars <- vapply(blocks, function(b) b$var, integer(1))
  tms <- .pr_funts_times(obj, idx)
  dn <- obj$dnames
  if (!is.null(dn) && is.list(dn)) dn <- dn[vars]
  vn <- obj$vnames
  if (!is.null(vn)) vn <- vn[vars]
  funts(X = Xl, basisobj = obj$basis[vars], argval = argv,
        method = if (grid_mode) "data" else "coefs",
        start = tms[1], end = tms[length(tms)], vnames = vn, dnames = dn, tname = obj$tname)
}

# convert a list of per-block d_b x n working series (centred or not; the
# caller adds the mean) at time indices idx to the class of the input
.pr_as_input <- function(dat, series, idx) {
  if (dat$type == "vector") {
    x <- as.numeric(series[[1]])
    tsp0 <- dat$meta$tsp
    if (!is.null(tsp0)) {
      x <- stats::ts(x, start = tsp0[1] + (idx[1] - 1) / tsp0[3], frequency = tsp0[3])
    }
    return(x)
  }
  if (dat$type == "matrix") {
    x <- do.call(cbind, lapply(series, function(s) as.numeric(s)))
    colnames(x) <- dat$meta$colnames
    tsp0 <- dat$meta$tsp
    if (!is.null(tsp0)) {
      x <- stats::ts(x, start = tsp0[1] + (idx[1] - 1) / tsp0[3], frequency = tsp0[3])
    }
    return(x)
  }
  els <- dat$meta$elements
  out <- lapply(els, function(el) {
    if (el$kind == "scalar") return(as.numeric(series[[el$blocks]]))
    .pr_make_funts(el$obj, dat$blocks[el$blocks], series[el$blocks], idx)
  })
  if (dat$type == "funts") return(out[[1]])
  names(out) <- vapply(els, function(el) if (is.null(el$name)) "" else el$name, character(1))
  out
}

# evaluation map from working coordinates (time vector of length d_tot) to the
# reporting points: identity for scalar and grid blocks, Phi_b R_b^{-1} for
# basis blocks (Phi_b = basis evaluated on the funts grid)
.pr_eval_map <- function(dat) {
  maps <- lapply(dat$blocks, function(blk) {
    if (blk$kind == "scalar" || isTRUE(blk$grid)) return(diag(blk$d))
    obj <- dat$meta$elements[[blk$element]]$obj
    obj$B_mat[[blk$var]] %*% blk$Rinv
  })
  nr <- vapply(maps, nrow, integer(1))
  nc <- vapply(maps, ncol, integer(1))
  E <- matrix(0, sum(nr), sum(nc))
  ro <- 0L
  co <- 0L
  for (b in seq_along(maps)) {
    E[ro + seq_len(nr[b]), co + seq_len(nc[b])] <- maps[[b]]
    ro <- ro + nr[b]
    co <- co + nc[b]
  }
  list(E = E, rows = split(seq_len(sum(nr)), rep(seq_along(nr), nr)))
}
