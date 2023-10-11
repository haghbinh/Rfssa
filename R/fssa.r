# --------------------------------------------------------------
#' Functional Singular Spectrum Analysis (FSSA)
#'
#' This function performs the decomposition (embedding and functional SVD steps)
#' for univariate (ufssa) or multivariate (mfssa) functional singular spectrum
#' analysis based on the input data type. The input can be a univariate or
#' multivariate functional time series (\code{\link{funts}}) object.
#'
#' @param Y an object of class \code{\link{funts}}.
#' @param L a positive integer, the window length.
#' @param ntriples a positive integer, the number of eigentriples for the
#' decomposition.
#' @param type a string indicating the type of FSSA: "ufssa" (default for
#' univariate FTS) or "mfssa" (default for multivariate FTS).
#'
#' @return An object of class \code{fssa}, containing functional objects,
#' eigenvalues, window length, and original data.
#'
#' @examples
#' \dontrun{
#' data("Callcenter")
#'
#' # FSSA Decomposition step:
#' L <- 28
#' U <- fssa(Callcenter, L)
#' plot(U, type = "values", d = 10)
#' plot(U, type = "vectors", d = 4)
#' plot(U, type = "paired", d = 6)
#' plot(U, type = "lcurves", d = 4, vars = 1)
#' plot(U, type = "lheats", d = 4)
#' plot(U, type = "wcor", d = 10)
#' plotly_funts(U$Lsingf[[1]])
#' plot(U$Lsingf[[2]])
#'
#' #--------------- Multivariate FSSA Example on bivariate -----------------------------
#' ## temperature curves and smoothed images of vegetation
#' data("Montana")
#'
#' # MFSSA Decomposition step:
#' L <- 45
#' U <- fssa(Montana, L)
#' plot(U, type = "values", d = 10)
#' plot(U, type = "vectors", d = 4)
#' plot(U, type = "lheats", d = 4)
#' plot(U, type = "lcurves", d = 4, vars = 1)
#' plot(U, type = "paired", d = 6)
#' plot(U, type = "periodogram", d = 4)
#' plot(U, type = "wcor", d = 10)
#' plotly_funts(U$Lsingf[[1]])
#' plot(U$Lsingf[[2]])
#'
#' }
#' @useDynLib Rfssa
#' @export
fssa <- function(Y, L = NA, ntriples = 20, type = "ufssa") {
  N <- Y$N
  if (is.na(L)) L <- floor(N / 2L)
  if (ntriples > L) {
    ntriples <- L
    warning("\"ntriples\" must be less than or equal to \"L\". Setting \"ntriples\" = \"L\"")
  }
  cat("Running, please wait...\n")
  p <- length(Y$dimSupp)
  if (p == 1 && type == "ufssa") {
    out <- ufssa(Y, L, ntriples)
  } else if (p > 1 || type == "mfssa") {
    out <- mfssa(Y, L, ntriples)
  } else {
    stop("Error in type or dimension.")
  }
  cat("Done.\n")
  class(out) <- "fssa"
  return(out)
}




#---------------------------------------------ufssa--------------------------------------------

# Embedding and decomposition stages of univariate functional singular spectrum analysis
ufssa <- function(Y, L, ntriples) {
  dimSupp <- Y$dimSupp
  N <- Y$N
  basis <- Y$B_mat[[1]]
  d <- ncol(Y$B_mat[[1]])
  grid <- as.matrix(Y$argval[[1]])
  K <- N - L + 1L
  if (dimSupp[[1]] == 1) {
    C_tilde <- t(onedG(A = basis %*% Y$coefs[[1]], B = basis, grid = grid))
    G <- onedG(A = basis, B = basis, grid = grid)
  } else {
    C_tilde <- t(twodG(A = basis %*% Y$coefs[[1]], B = basis, grid = grid))
    G <- twodG(A = basis, B = basis, grid = grid)
  }
  # Calculating Variance/Covariance Structure
  S0 <- SS(K, L, C_tilde, d)
  # Calculating Gram Matrix
  H <- CalculateInverse(Gram(K, L, G, d))
  # Calculating Eigen Triples
  Q <- eigs(AtimesB(H, S0), ntriples)
  # Returning results
  Q$values <- Re(Q$values)
  Q$vectors <- Re(Q$vectors)
  out <- list()
  Lsingf <- list()
  for (i in 1L:ntriples) {
    out[[i]] <- Y$B_mat[[1]] %*% Cofmat(d, L, Q$vectors[, i])
    if (Y$dimSupp[[1]] == 1) {
      X_mat <- out[[i]]
    } else {
      n1 <- length(unique(Y$argval[[1]][, 1]))
      n2 <- length(unique(Y$argval[[1]][, 2]))
      X_mat <- array(out[[i]], dim = c(n1, n2, L))
    }
    Lsingf[[i]] <- funts(X = X_mat, basisobj = Y$basis[[1]],
                         start = Y$time[1], end = Y$time[N], vnames = Y$vnames,
                         dnames = Y$dnames, tname = Y$tname)
  }
  out$values <- Q$values[1L:ntriples]
  out$L <- L
  out$N <- N
  out$Y <- Y
  out$RVectrs <- uV(out, ntriples)
  out$Lsingf <- Lsingf
  return(out)
}





#------------------------------------------------mfssa-----------------------------------------

# Embedding and decomposition stages of multivariate functional singular spectrum analysis.
mfssa <- function(Y, L, ntriples) {
  # get c plus plus code
  p <- length(Y$dimSupp)
  N <- Y$N
  C_tilde <- list()
  G_1 <- list()
  shifter <- matrix(nrow = 2, ncol = (p + 1L), data = 0L)
  Y_d <- matrix(data = 0, nrow = 1, ncol = p)
  # get inner product matrices
  for (j in 1:p) {
    grid <- as.matrix(Y$argval[[j]])
    if (Y$dimSupp[[j]] == 1) {
      C_tilde[[j]] <- t(onedG(A = Y$B_mat[[j]] %*% Y$coefs[[j]], B = Y$B_mat[[j]], grid = grid)) # old B
      G_1[[j]] <- t(onedG(A = Y$B_mat[[j]], B = Y$B_mat[[j]], grid = grid)) # old A
    } else {
      C_tilde[[j]] <- t(twodG(A = Y$B_mat[[j]] %*% Y$coefs[[j]], B = Y$B_mat[[j]], grid = grid))
      G_1[[j]] <- t(twodG(A = Y$B_mat[[j]], B = Y$B_mat[[j]], grid = grid))
    }
    shifter[1L, j + 1L] <- shifter[2L, j] + 1L
    shifter[2L, j + 1L] <- shifter[2L, j] + L * ncol(Y$B_mat[[j]])
    Y_d[1, j] <- ncol(Y$B_mat[[j]])
  }
  d <- cbind(0, Y_d * L)
  # Find the proper inner product matrices for j_k variables
  d_tilde <- sum(d) / L
  K <- N - L + 1L
  # Calculating Variance/Covariance Structure
  S0 <- SSM(K, L, d_tilde, p, C_tilde, shifter)
  # Calculating Gram Matrix
  H <- CalculateInverse(Gramm(K, L, p, d_tilde, G_1, shifter, d))
  # Calculating Eigen Triples
  Q <- eigs(AtimesB(H, S0), ntriples)
  # Returning results
  Q$values <- Re(Q$values)
  Q$vectors <- Re(Q$vectors)
  coefs0 <- Q$vectors
  p_c <- list()
  values <- Q$values[1L:ntriples]
  out <- list()
  Lsingf <- list()
  for (i in 1L:(ntriples)) {
    my_pcs <- list()
    X_mat <- list()
    for (j in 1L:p) {
      my_pcs[[j]] <- Y$B_mat[[j]] %*% Cofmat((d[j + 1L] / L), L, coefs0[(shifter[1L, (j + 1L)]:shifter[2L, (j + 1L)]), i])
      if (Y$dimSupp[[j]] == 1) {
        X_mat[[j]] <- my_pcs[[j]]
      } else {
        n1 <- length(unique(Y$argval[[j]][, 1]))
        n2 <- length(unique(Y$argval[[j]][, 2]))
        X_mat[[j]] <- array(my_pcs[[j]], dim = c(n1, n2, L))
      }
    }
    Lsingf[[i]] <- funts(X = X_mat, basisobj = Y$basis,
                         start = Y$time[1], end = Y$time[N], vnames = Y$vnames,
                         dnames = Y$dnames, tname = Y$tname)
    out[[i]] <- my_pcs
  }
  out$values <- values
  out$L <- L
  out$N <- N
  out$Y <- Y
  out$RVectrs <- mV(out, ntriples)
  out$Lsingf <- Lsingf
  return(out)
}
