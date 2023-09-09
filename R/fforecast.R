#--------------------------------------------------------------
#' Functional Singular Spectrum Analysis Recurrent Forecasting and Vector Forecasting
#'
#' This function performs functional singular spectrum analysis (FSSA) recurrent forecasting (FSSA R-forecasting) or vector forecasting (FSSA V-forecasting) of univariate or multivariate functional time series (\code{\link{fts}}) observed over a one-dimensional domain.
#' @return A list of objects of class \code{\link{fts}} where each fts corresponds to a forecasted group.
#' @param U An object of class \code{\link{fssa}} that holds the decomposition.
#' @param groups A list of numeric vectors where each vector includes indices of elementary components of a group used for reconstruction and forecasting.
#' @param h An integer that specifies the forecast horizon.
#' @param method A character string specifying the type of forecasting to perform either \code{"recurrent"} or \code{"vector"}.
#' @param tol A double specifying the amount of tolerated error in the approximation of the matrix that corresponds with the operator formed using a Neumann series leveraged in both forecasting algorithms.
#' @examples
#' \dontrun{
#' data("Callcenter")
#' L <- 28
#' U <- fssa(Callcenter, L)
#' groups <- list(1,1:7)
#' ## Perform FSSA R-forecast
#' pr_R <- fforecast(U = U, groups = groups, h = 30, method = "recurrent")
#'
#' plotly_funts(pr_R[[1]], main = "Call Center Mean Component Recurrent Forecast",
#'              xlab = "Time (6 minutes aggregated)",
#'              ylab = "Sqrt of Call Numbers",type="line",
#'              xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),xticklocs =
#'                list(c(1,60,120,180,240)))
#' plotly_funts(pr_R[[2]], main = "Call Center Recurrent Forecast from 7'th first components",
#'              xlab = "Time (6 minutes aggregated)",
#'              ylab = "Sqrt of Call Numbers",type="line",
#'              xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),xticklocs =
#'                list(c(1,60,120,180,240)))
#'
#'
#' ## Perform FSSA V-forecast
#' pr_V <- fforecast(U = U, groups = groups, h = 30, method = "vector", tol = 10^-3)
#'
#' plotly_funts(pr_V[[1]], mains = "Call Center Mean Component Vector Forecast",
#'      xlabels = "Time (6 minutes aggregated)",
#'      ylabels = "Sqrt of Call Numbers",type="line",
#'      xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),xticklocs =
#'        list(c(1,60,120,180,240)))
#' plotly_funts(pr_V[[2]], mains = "Call Center Vector Forecast from 7'th first components",
#'      xlabels = "Time (6 minutes aggregated)",
#'      ylabels = "Sqrt of Call Numbers",type="line",
#'      xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),xticklocs =
#'        list(c(1,60,120,180,240)))
#' }
#'
#' @export
fforecast <- function(U, groups = list(1), h = 1, method = "recurrent", tol = 10^-3) {
  for (j in 1:length(U$Y$coefs)) {
    if (U$Y$dimSupp[[j]] > 1) {
      stop("Current forecasting routines only support fts whose variables are observed over one-dimensional domains. Forecasting of fts variables whose domains have dimension greater than one is under development.")
    }
  }
  cat("Running, please wait...\n")
  if (class(U[[1]])[[1]] != "list") {
    out <- ufforecast(U = U, groups = groups, h = h, method = method, tol = tol)
  } else {
    out <- mfforecast(U = U, groups = groups, h = h, method = method, tol = tol)
  }
  cat("Done.\n")
  return(out)
}


#------------------------------ufforecast-----------------------------------------------------------------
# FSSA Recurrent and Vector Forecasting of univariate FTS

ufforecast <- function(U, groups = list(c(1)), h = 1, method = "recurrent", tol = 10^-3) {
  out <- list()
  Y <- U$Y
  basisobj <- Y$basis[[1]]
  time_st <- Y$time[1]
  time_en <- Y$time[N]
  for (a in 1:length(groups)) {
    g <- groups[[a]]
    # Define prediction space
    basis <- U$Y$B_mat[[1]]
    N <- U$N
    d <- ncol(basis)
    L <- U$L
    K <- N - L + 1
    k <- length(g)
    G <- t(basis) %*% basis
    G_inv <- solve(G)
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
      fssa_fore <- matrix(data = 0, nrow = d, ncol = h)
      for (m in 1:h) {
        for (j in 1:(L - 1)) {
          E_j <- matrix(data = NA, nrow = d, ncol = k)
          for (n in 1:k) {
            E_j[, n] <- G_inv %*% t(basis) %*% U[[g[n]]][, j]
          }
          A_j <- Neu %*% D %*% t(E_j)
          fssa_fore[, m] <- fssa_fore[, m] + A_j %*% as.matrix(Q[[1]]$coefs[[1]][, (N + j - L + m)])
        }
        Q[[1]]$coefs[[1]] <- cbind(Q[[1]]$coefs[[1]], fssa_fore[, m])
      }
      funts_out <- funts(X = (basis %*% fssa_fore), basisobj = basisobj, argval = Q[[1]]$argval[[1]], start = time_st, end = time_en)
      out[[a]] <- funts_out
    } else if (method == "vector") {
      # FSSA V-forecasting
      F_matrix <- matrix(data = NA, nrow = ((L - 1) * d), ncol = k)
      for (n in 1:k) {
        F_matrix[, n] <- matrix(data = G_inv %*% t(basis) %*% U[[g[n]]][, 1:(L - 1)], nrow = ((L - 1) * d), ncol = 1)
      }
      P <- F_matrix %*% t(F_matrix) + F_matrix %*% t(D) %*% Neu %*% D %*% t(F_matrix)
      S <- array(data = 0, dim = c(d, (K + h), L))
      for (j in 1L:length(g)) {
        S[, (1:K), ] <- S[, (1:K), ] + ufproj(U, g[j], d)
      }
      for (m in 1:h) {
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
      predictions <- S[, (K + 1):(K + h), L]
      funts_out <- funts(X = (basis %*% predictions), basisobj = basisobj, argval = U$Y$argval[[1]], start = time_st, end = time_en)
      out[[a]] <- funts_out
    }
  }
  return(out)
}


#------------------------------mfforecast-----------------------------------------------------------------

# FSSA Recurrent and Vector Forecasting of Multivariate FTS

mfforecast <- function(U, groups = list(c(1)), h = 1, method = "recurrent", tol = 10^-3) {
  out <- list()
  for (a in 1:length(groups)) {
    g <- groups[[a]]
    # Define prediction space
    basis <- list()
    G_inv <- list()
    p <- length(U$Y$coefs)
    d <- matrix(data = 0, nrow = 1, ncol = p)
    N <- U$N
    L <- U$L
    K <- N - L + 1
    shifter <- matrix(data = 0, nrow = 2, ncol = (p+1))
    k <- length(g)
    D <- matrix(data = NA, nrow = 1, ncol = k)
    for (j in 1:p) {
      basis[[j]] <- U$Y$B_mat[[j]]
      G_inv[[j]] <- solve(t(basis[[j]])%*%basis[[j]])
      d[j] <- ncol(U$Y$B_mat[[j]])
      shifter[1, j+1] <- shifter[2, j] + 1
      shifter[2, j+1] <- shifter[2, j] + ncol(U$Y$B_mat[[j]])
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
      fssa_fore <- matrix(data = 0, nrow = sum(d), ncol = h)
      for (m in 1:h) {
        for (j in 1:(L - 1)) {
          my_obs <- matrix(data = 0, nrow = sum(d), ncol = 1)
          E_j <- matrix(data = NA, nrow = sum(d), ncol = k)
          for (q in 1:p) {
            for (n in 1:k) {
              E_j[(shifter[1, (q+1)]:shifter[2, (q+1)]), n] <- G_inv[[q]] %*% t(basis[[q]]) %*% U[[g[n]]][[q]][, j]
            }
            my_obs[(shifter[1, (q+1)]:shifter[2, (q+1)]), 1] <- Q[[1]]$coefs[[q]][, (N + j - L + m)]
          }
          A_j <- Neu %*% D %*% t(E_j)
          fssa_fore[, m] <- fssa_fore[, m] + A_j %*% my_obs
        }
        for(q in 1:p) Q[[1]]$coefs[[q]] <- cbind(Q[[1]]$coefs[[q]], fssa_fore[(shifter[1, (q+1)]:shifter[2, (q+1)]), m]);
      }
      for (q in 1:p) {
        out_g[[q]] <- basis[[q]] %*% fssa_fore[(shifter[1, (q+1)]:shifter[2, (q+1)]), ]
      }
      fts_out <- Rfssa::fts(out_g, basis, U$Y$argval)
      fts_out$B_matasis_type <- U$Y$B_matasis_type
      out[[a]] <- fts_out
    } else if (method == "vector") {
      # MFSSA V-forecasting
      out_g <- list()
      Y <- matrix(data = NA, nrow = ((L - 1) * sum(d)), ncol = k)
      Lshifter <- matrix(data = 0, nrow = 2, ncol = (p+1))
      for (j in 1:p) {
        Lshifter[1, (j+1)] <- Lshifter[2, j] + 1
        Lshifter[2, (j+1)] <- Lshifter[2, j] + (L-1)*ncol(U$Y$B_mat[[j]])
        for (n in 1:k) {
          Y[(Lshifter[1, (j+1)]:Lshifter[2, (j+1)]), n] <- matrix(data = G_inv[[j]] %*% t(basis[[j]]) %*% U[[g[n]]][[j]][, 1:(L - 1)], nrow = ((L - 1) * d[j]), ncol = 1)
        }
      }
      P <- Y %*% t(Y) + Y %*% t(D) %*% Neu %*% D %*% t(Y)
      S <- list()
      for(j in 1:p) S[[j]] <- array(data = 0, dim = c(d[j], (K + h), L));
      for (k in 1L:length(g)) {
        projection <- mfproj(U, g[k])
        for (j in 1:p) {
          S[[j]][, 1:K, ] <- S[[j]][, 1:K, ] + projection[[j]]
        }
      }
      for (m in 1:h) {
        obs <- matrix(data = NA, nrow = ((L - 1) * sum(d)), ncol = 1)
        for (j in 1:p) {
          obs[(Lshifter[1, (j+1)]:Lshifter[2, (j+1)]), 1] <- matrix(data = S[[j]][, (K + (m - 1)), 2:L], nrow = ((L - 1) * d[j]), ncol = 1)
        }
        pr <- P %*% obs
        pr_mat <- matrix(data = NA, nrow = sum(d), ncol = L - 1)
        pr_1 <- matrix(data = 0, nrow = sum(d), ncol = 1)
        for (j in 1:(L - 1)) {
          my_obs <- matrix(data = 0, nrow = sum(d), ncol = 1)
          E_j <- matrix(data = NA, nrow = sum(d), ncol = k)
          for (q in 1:p) {
            pr_mat[(shifter[1, (q+1)]:shifter[2, (q+1)]), ] <- pr[(Lshifter[1, (q+1)]:Lshifter[2, (q+1)]), 1]
            for (n in 1:k) {
              E_j[(shifter[1, (q+1)]:shifter[2, (q+1)]), n] <- G_inv[[q]] %*% t(basis[[q]]) %*% U[[g[n]]][[q]][, j]
            }
            my_obs <- pr_mat[, j]
          }
          A_j <- Neu %*% D %*% t(E_j)
          pr_1 <- pr_1 + A_j %*% my_obs
        }
        for(q in 1:p) S[[q]][, (K + m), ] <- matrix(data = cbind(pr_mat, pr_1)[(shifter[1, (q+1)]:shifter[2, (q+1)]), ], nrow = d[q], ncol = L);
      }

      for (q in 1:p) {
        S[[q]] <- fH(S[[q]], d[q])
        out_g[[q]] <- basis[[q]] %*% S[[q]][, (K + 1):(K + h), L]
      }

      fts_out <- Rfssa::fts(out_g, basis, U$Y$argval)
      fts_out$B_matasis_type <- U$Y$B_matasis_type
      out[[a]] <- fts_out
    }
  }

  return(out)
}
