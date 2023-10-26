#--------------------------------------------------------------
#' Functional Singular Spectrum Analysis Recurrent and Vector Forecasting
#'
#' Perform functional singular spectrum analysis (FSSA) recurrent forecasting
#' (FSSA R-forecasting) or vector forecasting (FSSA V-forecasting) on univariate
#' or multivariate functional time series (\code{\link{funts}}) observed over a
#' one-dimensional domain.
#'
#' @param U an object of class \code{\link{fssa}} holding the decomposition.
#' @param groups a list of numeric vectors where each vector includes indices
#'               of elementary components of a group used for reconstruction and
#'               forecasting.
#' @param len integer, the desired length of the forecasted FTS.
#' @param method a character string specifying the type of forecasting to perform:
#'              - "recurrent" for FSSA R-forecasting.
#'              - "vector" for FSSA V-forecasting.
#' @param only.new logical, if 'TRUE' then only forecasted FTS are returned, whole FTS otherwise.
#' @param tol a double specifying the tolerated error in the approximation of
#'            the matrix used in forecasting algorithms.
#'
#' @return An object of class `fforecast` which is a list of objects of class
#'         \code{\link{funts}}, where each one corresponds to a forecasted group.
#'
#' @examples
#' \dontrun{
#' data("Callcenter")
#' U <- fssa(Callcenter, L = 28)
#' groups <- list(1, 1:7)
#'
#' ## Perform FSSA R-forecast
#' pr_R <- fforecast(
#'   U = U, groups = groups, only.new = FALSE,
#'   len = 30, method = "recurrent"
#' )
#'
#' plot(pr_R,  group_index = 1 )
#'
#'
#' plotly_funts(pr_R[[2]], main = "group = '1:7'")
#'
#' ## Perform FSSA V-forecast
#' pr_V <- fforecast(U = U, groups = groups, len= 30, method = "vector")
#'
#' plot(pr_V, group_index = 1)
#'
#' plotly_funts(pr_V[[2]], type = "3Dline" , main = "group = '1:7'")
#'
#' # Multivariate forecasting example:
#' data("Montana")
#' time <- Montana$time
#' grid <- list(0:23, list(1:33, 1:33))
#' montana <- eval.funts(Montana, argvals = grid)
#' montana[[2]] <- array(
#'   scale(montana[[2]][, , ],
#'     center = min(montana[[2]][, , ]),
#'     scale = max(montana[[2]][, , ]) - min(montana[[2]][, , ])
#'   ),
#'   dim = c(33, 33, 133)
#' )
#' ## Kernel density estimation of pixel intensity
#' NDVI <- matrix(NA, nrow = 512, ncol = 133)
#' for (i in 1:133) NDVI[, i] <- (density(montana[[2]][, , i], from = 0, to = 1)$y)
#'
#' ## Define functional objects
#' bs1 <- Montana$basis[[1]]
#'
#' require(fda)
#' bs2 <- create.bspline.basis(nbasis = 15)
#' Y <- funts(X = list(montana[[1]], NDVI), basisobj = list(bs1, bs2),
#'             vnames = c("Temperature", "NDVI Density"),
#'             dnames = c("Time", "NDVI"),
#'             tname = "Date")
#'
#' plotly_funts(Y,
#'   main = c("Temperature", "NDVI"),
#'   xticklocs = list(c(0, 6, 12, 18, 23), seq(1, 512, len = 9)),
#'   xticklabels = list(c(0, 6, 12, 18, 23), seq(0, 1, len = 9))
#' )
#'
#' U <- fssa(Y = Y, L = 45)
#' plotly_funts(U$Lsingf[[1]])
#' plot(U$Lsingf[[2]])
#'
#' groups <- list(1, 1:3)
#' pr_R <- fforecast(U = U, groups = groups,
#'                    only.new = FALSE, len = 10, method = "recurrent")
#' plot(pr_R)
#' plotly_funts(pr_R[[2]], main = "Recurrent method, group = '1:3'")
#'
#' pr_V <- fforecast(U = U, groups = groups, len = 10, method = "vector")
#' plot(pr_V, group_index = 1)
#' plotly_funts(pr_V[[2]], main = "Vector method, group = '1:3'")
#' }
#'
#' @export
fforecast <- function(U, groups, len = 1, method = "recurrent", only.new = TRUE, tol = NULL) {
  if (is.null(tol)) tol <- 10^(-3)
  for (j in 1:length(U$Y$coefs)) {
    if (U$Y$dimSupp[[j]] > 1) {
      stop("Current forecasting routines only support fts whose variables are observed over one-dimensional domains. Forecasting of fts variables whose domains have dimension greater than one is under development.")
    }
  }
  if (is.numeric(groups) & !is.list(groups)) groups <- list(groups)
  cat("Running, please wait...\n")
  if (class(U[[1]])[[1]] != "list") {
    out <- ufforecast(U = U, groups = groups, len = len, method = method, only.new = only.new, tol = tol)
  } else {
    out <- mfforecast(U = U, groups = groups, len = len, method = method, only.new = only.new, tol = tol)
  }
  cat("Done.\n")
  Y <- U$Y
  N <- U$N
  time_dif <- Y$time[N] - Y$time[N-1]
  time_st <- (Y$time[N] + time_dif)
  time_en <- Y$time[N] + len * time_dif
  out$predicted_time <- seq(from = time_st, to = time_en, length.out = len)
  out$groups <- groups
  out$original_funts <- U$Y
  out$method <- method
  class(out) <- "fforecast"
  return(out)
}


#------------------------------ufforecast-----------------------------------------------------------------
# FSSA Recurrent and Vector Forecasting of univariate FTS

ufforecast <- function(U, groups, len = 1, method = "recurrent", only.new = TRUE, tol = 10^-3) {
  out <- list()
  Y <- U$Y
  N <- Y$N
  basis <- U$Y$B_mat[[1]]
  d <- ncol(basis)
  L <- U$L
  K <- N - L + 1
  G <- t(basis) %*% basis
  G_inv <- solve(G)
  basisobj <- Y$basis[[1]]
  time_dif <- Y$time[N] - Y$time[N-1]
  if (only.new) time_st <- (Y$time[N] + time_dif) else time_st <- Y$time[1]
  time_en <- Y$time[N] + len * time_dif
  for (a in 1:length(groups)) {
    g <- groups[[a]]
    # Define prediction space
    k <- length(g)
    D <- matrix(data = NA, nrow = d, ncol = k)
    for (n in 1:k) D[, n] <- G_inv %*% t(basis) %*% U[[g[n]]][, L]

    # Define Truncated Neumann Series
    NU_1 <- D %*% t(D)
    NU <- list()
    NU[[1]] <- diag(1, d)
    norm <- 1
    l <- 0
    while (norm > tol) {
      l <- l + 1
      norm <- norm(NU[[l]])
      NU[[(l + 1)]] <- NU_1 %*% NU[[l]]
    }
    Neu <- NU[[1]]
    for (n in 2:l) Neu <- Neu + NU[[n]]

    if (method == "recurrent") {
      # FSSA R-forecasting
      # Reconstruct signal
      Q <- freconstruct(U, groups = list(g))
      fore_r <- matrix(data = 0, nrow = d, ncol = len )
      for (m in 1:len ) {
        for (j in 1:(L - 1)) {
          E_j <- matrix(data = NA, nrow = d, ncol = k)
          for (n in 1:k) {
            E_j[, n] <- G_inv %*% t(basis) %*% U[[g[n]]][, j]
          }
          A_j <- Neu %*% D %*% t(E_j)
          fore_r[, m] <- fore_r[, m] + A_j %*% as.matrix(Q[[1]]$coefs[[1]][, (N + j - L + m)])
        }
        Q[[1]]$coefs[[1]] <- cbind(Q[[1]]$coefs[[1]], fore_r[, m])
      }
      if (!only.new) {
        fore_r <- cbind(U$Y$coefs[[1]], fore_r)
      }
      funts_out <- funts(X = (basis %*% fore_r), basisobj = basisobj, argval = Y$argval[[1]],
                         start = time_st, end = time_en, vnames = Y$vnames, dnames = Y$dnames, tname = Y$tname)
      out[[a]] <- funts_out
    } else if (method == "vector") {
      # FSSA V-forecasting
      F_matrix <- matrix(data = NA, nrow = ((L - 1) * d), ncol = k)
      for (n in 1:k) {
        F_matrix[, n] <- matrix(data = G_inv %*% t(basis) %*% U[[g[n]]][, 1:(L - 1)], nrow = ((L - 1) * d), ncol = 1)
      }
      P <- F_matrix %*% t(F_matrix) + F_matrix %*% t(D) %*% Neu %*% D %*% t(F_matrix)
      S <- array(data = 0, dim = c(d, (K + len ), L))
      for (j in 1L:length(g)) {
        S[, (1:K), ] <- S[, (1:K), ] + ufproj(U, g[j], d)
      }
      for (m in 1:len) {
        obs <- matrix(data = S[, (K + (m - 1)), 2:L], nrow = ((L - 1) * d), ncol = 1)
        pr <- P %*% obs
        pr <- matrix(data = pr, nrow = d, ncol = (L - 1))
        pr_1 <- matrix(data = 0, nrow = d, ncol = 1)
        for (j in 1:(L - 1)) {
          E_j <- matrix(data = NA, nrow = d, ncol = k)
          for (n in 1:k) {
            E_j[, n] <- G_inv %*% t(basis) %*% U[[g[n]]][, j]
          }
          A_j <- Neu %*% D %*% t(E_j)

          pr_1[, 1] <- pr_1[, 1] + A_j %*% S[, (K + m - 1), (j + 1)]
        }
        pr <- cbind(pr, pr_1)
        S[, (K + m), ] <- pr
      }

      S <- fH(S, d)
      fore_v <- S[, (K + 1):(K + len), L]
      if (!only.new) {
        fore_v <- cbind(U$Y$coefs[[1]], fore_v)
      }
      funts_out <- funts(X = (basis %*% fore_v), basisobj = basisobj, argval = Y$argval[[1]],
                         start = time_st, end = time_en, vnames = Y$vnames, dnames = Y$dnames, tname = Y$tname)
      out[[a]] <- funts_out
    }
  }
  return(out)
}


#------------------------------mfforecast-----------------------------------------------------------------

# FSSA Recurrent and Vector Forecasting of Multivariate FTS

mfforecast <- function(U, groups, len = 1, method = "recurrent", only.new = TRUE, tol = 10^-3) {
  out <- list()
  Y <- U$Y
  N <- Y$N
  basis <- U$Y$B_mat
  G_inv <- lapply(basis, function(x) {
    solve(t(x) %*% x)
  })
  d <- sapply(basis, ncol)
  p <- length(Y$dimSupp)
  L <- U$L
  K <- N - L + 1
  shifter <- matrix(data = 0, nrow = 2, ncol = (p + 1))
  for (j in 1:p) {
    shifter[1, j + 1] <- shifter[2, j] + 1
    shifter[2, j + 1] <- shifter[2, j] + ncol(U$Y$B_mat[[j]])
  }
  basisobj <- Y$basis
  time_dif <- Y$time[N] - Y$time[N-1]
  if (only.new) time_st <- (Y$time[N] + time_dif) else time_st <- Y$time[1]
  time_en <- Y$time[N] + len * time_dif

  for (a in 1:length(groups)) {
    g <- groups[[a]]
    k <- length(g)
    D <- matrix(data = NA, nrow = 1, ncol = k)
    for (j in 1:p) {
      D_j <- matrix(data = NA, nrow = d[j], ncol = k)
      for (n in 1:k) D_j[, n] <- solve(t(basis[[j]]) %*% basis[[j]]) %*% t(basis[[j]]) %*% U[[g[n]]][[j]][, L]
      D <- rbind(D, D_j)
    }
    D <- D[(2:nrow(D)), ]

    # Define Truncated Neumann Series
    NU_1 <- D %*% t(D)
    NU <- list()
    NU[[1]] <- diag(1, sum(d))
    norm <- 1
    l <- 0
    while (norm > tol) {
      l <- l + 1
      norm <- norm(NU[[l]])
      NU[[(l + 1)]] <- NU_1 %*% NU[[l]]
    }
    Neu <- NU[[1]]
    for (n in 2:l) Neu <- Neu + NU[[n]]

    if (method == "recurrent") {
      # MFSSA R-forecasting
      out_g <- list()
      Q <- freconstruct(U, groups = list(g))
      fore_r <- matrix(data = 0, nrow = sum(d), ncol = len)
      for (m in 1:len) {
        for (j in 1:(L - 1)) {
          my_obs <- matrix(data = 0, nrow = sum(d), ncol = 1)
          E_j <- matrix(data = NA, nrow = sum(d), ncol = k)
          for (q in 1:p) {
            for (n in 1:k) {
              E_j[(shifter[1, (q + 1)]:shifter[2, (q + 1)]), n] <- G_inv[[q]] %*% t(basis[[q]]) %*% U[[g[n]]][[q]][, j]
            }
            my_obs[(shifter[1, (q + 1)]:shifter[2, (q + 1)]), 1] <- Q[[1]]$coefs[[q]][, (N + j - L + m)]
          }
          A_j <- Neu %*% D %*% t(E_j)
          fore_r[, m] <- fore_r[, m] + A_j %*% my_obs
        }
        for (q in 1:p) Q[[1]]$coefs[[q]] <- cbind(Q[[1]]$coefs[[q]], fore_r[(shifter[1, (q + 1)]:shifter[2, (q + 1)]), m])
      }
      for (q in 1:p) {
        if (only.new){
          X0 <- fore_r[(shifter[1, (q + 1)]:shifter[2, (q + 1)]), ]
        } else {
          X0 <- cbind(U$Y$coefs[[q]], fore_r[(shifter[1, (q + 1)]:shifter[2, (q + 1)]), ])
        }
        out_g[[q]] <- basis[[q]] %*% X0
      }
      funts_out <- funts(X = out_g, basisobj = basisobj, argval = Y$argval,
                         start = time_st, end = time_en, vnames = U$Y$vnames,
                         dnames = U$Y$dnames, tname = U$Y$tname)
      out[[a]] <- funts_out
    } else if (method == "vector") {
      # MFSSA V-forecasting
      out_g <- list()
      Y <- matrix(data = NA, nrow = ((L - 1) * sum(d)), ncol = k)
      Lshifter <- matrix(data = 0, nrow = 2, ncol = (p + 1))
      for (j in 1:p) {
        Lshifter[1, (j + 1)] <- Lshifter[2, j] + 1
        Lshifter[2, (j + 1)] <- Lshifter[2, j] + (L - 1) * ncol(U$Y$B_mat[[j]])
        for (n in 1:k) {
          Y[(Lshifter[1, (j + 1)]:Lshifter[2, (j + 1)]), n] <- matrix(data = G_inv[[j]] %*% t(basis[[j]]) %*% U[[g[n]]][[j]][, 1:(L - 1)], nrow = ((L - 1) * d[j]), ncol = 1)
        }
      }
      P <- Y %*% t(Y) + Y %*% t(D) %*% Neu %*% D %*% t(Y)
      S <- list()
      for (j in 1:p) S[[j]] <- array(data = 0, dim = c(d[j], (K + len), L))
      for (k in 1L:length(g)) {
        projection <- mfproj(U, g[k])
        for (j in 1:p) {
          S[[j]][, 1:K, ] <- S[[j]][, 1:K, ] + projection[[j]]
        }
      }
      for (m in 1:len) {
        obs <- matrix(data = NA, nrow = ((L - 1) * sum(d)), ncol = 1)
        for (j in 1:p) {
          obs[(Lshifter[1, (j + 1)]:Lshifter[2, (j + 1)]), 1] <- matrix(data = S[[j]][, (K + (m - 1)), 2:L], nrow = ((L - 1) * d[j]), ncol = 1)
        }
        pr <- P %*% obs
        pr_mat <- matrix(data = NA, nrow = sum(d), ncol = L - 1)
        pr_1 <- matrix(data = 0, nrow = sum(d), ncol = 1)
        for (j in 1:(L - 1)) {
          my_obs <- matrix(data = 0, nrow = sum(d), ncol = 1)
          E_j <- matrix(data = NA, nrow = sum(d), ncol = k)
          for (q in 1:p) {
            pr_mat[(shifter[1, (q + 1)]:shifter[2, (q + 1)]), ] <- pr[(Lshifter[1, (q + 1)]:Lshifter[2, (q + 1)]), 1]
            for (n in 1:k) {
              E_j[(shifter[1, (q + 1)]:shifter[2, (q + 1)]), n] <- G_inv[[q]] %*% t(basis[[q]]) %*% U[[g[n]]][[q]][, j]
            }
            my_obs <- pr_mat[, j]
          }
          A_j <- Neu %*% D %*% t(E_j)
          pr_1 <- pr_1 + A_j %*% my_obs
        }
        for (q in 1:p) S[[q]][, (K + m), ] <- matrix(data = cbind(pr_mat, pr_1)[(shifter[1, (q + 1)]:shifter[2, (q + 1)]), ], nrow = d[q], ncol = L)
      }

      for (q in 1:p) {
        S[[q]] <- fH(S[[q]], d[q])

        if (only.new){
          X0 <- S[[q]][, (K + 1):(K + len), L]
        } else {
          X0 <- cbind(U$Y$coefs[[q]], S[[q]][, (K + 1):(K + len), L])
        }
        out_g[[q]] <- basis[[q]] %*% X0
      }
      funts_out <- funts(X = out_g, basisobj = basisobj, argval = U$Y$argval,
                         start = time_st, end = time_en, vnames = U$Y$vnames,
                         dnames = U$Y$dnames, tname = U$Y$tname)
      out[[a]] <- funts_out
    }
  }
  return(out)
}


# =======================================================================
#'
#' Plot Method for FSSA Forecast (fforecast) Class
#'
#' Create visualizations of FSSA Forecast (fforecast) class. This function supports
#' plotting `fforecast` data with one-dimensional or two-dimensional domains.
#'
#' @param x an object of class \code{fforecast}.
#' @param group_index an integer specifying the group index for the plot.
#' @param ask logical: If `TRUE`, and  `group_index` be `NULL`, after printing the first grouping graphic, it will pause when the user asks for the next group graphic and wait.
#' @param npts number of grid points for the plots.
#' @param obs observation number (for two-dimensional domains).
#' @param main main title for the plot.
#' @param col specify the predicted FTS color; if it is `NULL`, it will be set as the default.
#' @param ori_col specify the original FTS color; if it is `NULL`, it will be set as the default.
#' @param type type of plot ("l" for line, "p" for points, etc.).
#' @param lty line type (1 for solid, 2 for dashed, etc.).
#' @param ... additional graphical parameters passed to plotting functions.
#' @seealso \code{\link{fforecast}}
#'
#' @examples
#' \dontrun{
#' # Example with one-dimensional domain
#' data("Callcenter")
#' # FSSA Decomposition step:
#' fssa_results <- fssa(Callcenter, L = 28)
#'
#' # Perform FSSA R-forecasting
#' pr_V <- fforecast(U = fssa_results, groups = list(1,1:7),
#'                   len = 14, method = "vector", only.new = FALSE)
#'
#' plot(pr_V)
#' }
#'
#' @export
plot.fforecast <- function(x, group_index = NULL, ask = TRUE, npts = 100, obs = 1,
                           main = NULL, col = NULL, ori_col = NULL, type = "l",
                           lty = 1, ...) {
  old <- par()
  exclude_pars <- c("cin", "cra", "csi", "cxy", "din", "page")
  ind <- which(!(names(old) %in% exclude_pars))
  on.exit(par(old[ind]))
  N <- x$original_funts$N
  h <- length(x$predicted_time)
  if(is.null(ori_col)) ori_col <- rep('snow3', N)
  if(is.null(col)) col <- rep("deepskyblue4", h)
  N1 <- x[[1]]$N
  if (N1 == h) col1 <- col else  col1 <- c(ori_col, col)
  flag <- FALSE
  if (is.null(group_index)) {
    group_index <- 1:length(x$groups)
    if (length(x$groups) > 1) flag <- TRUE
  }
  if (flag) {
    par(ask = ask)
    for (ipc in group_index) {
      obj <- x[[ipc]]
      plot(obj,col = col1, npts = npts, obs = obs, main = paste(main, "Group index:", ipc), type = type, lty = lty, ...)
    }
  } else {
    obj <- x[[group_index]]
    plot(obj,col = col1, npts = npts, obs = obs, main = paste(main, "Group index:", group_index), type = type, lty = lty, ...)
  }
}



# =======================================================================
#' Custom Print Method for FSSA Forecast (fforecast) class
#'
#' This custom print method is designed for objects of the FSSA Forecast (fforecast) class.
#' It provides a summary of the fforecast object.
#'
#' @param x an object of class "fforecast" to be printed.
#' @param ...	 further arguments passed to or from other methods.
#'
#' @examples
#' \dontrun{
#' # Example with one-dimensional domain
#' data("Callcenter")
#' # FSSA Decomposition step:
#' fssa_results <- fssa(Callcenter, L = 28)
#'
#' # Perform FSSA R-forecasting
#' pr_R <- fforecast(U = fssa_results,
#'                   groups = c(1:3),
#'                   len = 14,
#'                   method = "recurrent")
#' print(pr_R)
#'}
#' @export
print.fforecast <- function(x, ...) {
  cat("\nFSSA Forecast (fforecast) class:")
  cat("\nGroups: ")
  cat(str(x$groups))
  cat("Prediction method: ", x$method)
  cat("\nPredicted series length: ", length(x$predicted_time))
  cat("\nPredicted time: ")
  cat(str(x$predicted_time))
  cat("\n---------The original series-----------")
  print(x$original_funts)
}
