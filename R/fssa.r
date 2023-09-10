#' Functional Singular Spectrum Analysis
#'
#' This is a function which performs the decomposition (including embedding
#'  and  functional SVD steps) stage for univariate functional singular spectrum analysis (ufssa)
#'  or multivariate functional singular spectrum analysis (mfssa). The algorithm (ufssa or mfssa) is chosen based on
#'  whether the supplied input is a univariate or
#'  multivariate functional time series (\code{\link{fts}}) object. The \code{type} parameter can also be set to \code{"mfssa"} if the user wishes to perform ufssa of a univariate \code{\link{fts}} object using the mfssa code. Also note that the variables of the \code{\link{fts}} maybe observed over different dimensional domains where the maximum dimension currently supported is two.
#' @return An object of class \code{fssa}, which is a list of functional objects and the following components:
#' \item{values}{A numeric vector of eigenvalues.}
#' \item{L}{The specified window length.}
#' \item{N}{The length of the functional time series.}
#' \item{Y}{The original functional time series.}
#' @param Y An object of class \code{\link{fts}}.
#' @param L A positive integer giving the window length.
#' @param ntriples A positive integer specifying the number of eigentriples to calculate in the decomposition.
#' @param type A string indicating which type of fssa to perform. Use \code{type="ufssa"} to perform univariate fssa (default for univariate fts). Use \code{type="mfssa"} to perform multivariate fssa (default for multivariate fts).
#' @importFrom RSpectra eigs
#' @examples
#' \dontrun{
#' #---------------- Univariate FSSA Example on Callcenter data-------------------
#' data("Callcenter")
#' plot(Callcenter,lwd=2, npts = 200, col = "deepskyblue4",
#'      main = "Call Center Data",
#'      xlab = "Time (6 minutes aggregated)",
#'      ylab = "Sqrt of Call Numbers")
#'
#' plotly_funts(Callcenter,
#'              main = "Call Center Data Line Plot",
#'              xlab = "Time (6 minutes aggregated)",
#'              ylab = "Sqrt of Call Numbers",type="line",
#'             xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),
#'              xticklocs = list(c(1,60,120,180,240)))
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
#'
#' plotly_funts(U$Lsingf[[1]])
#' plot(U$Lsingf[[2]])
#'
#' # FSSA Reconstruction step:
#' gr <- list(1, 2:3, 4:5, 6:7, 1:7)
#' Q <- freconstruct(U, gr)
#' plotly_funts(Q[[1]], mains = "Call Center Mean Component",
#'              xlab = "Time (6 minutes aggregated)",
#'              ylab = "Sqrt of Call Numbers",type="line",
#'              xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),xticklocs =
#'                list(c(1,60,120,180,240)))
#' plot(Q[[2]], main = "Call Center First Periodic Component",
#'      xlab = "Time (6 minutes aggregated)",
#'      ylab = "Sqrt of Call Numbers")
#'
#' #--------------- Multivariate FSSA Example on bivariate -----------------------------
#' ## temperature curves and smoothed images of vegetation
#' data("Montana")
#' plot(Montana, obs = 2,
#'      main = c("Temperature Curves", "NDVI Images,"),
#'      xlab = c("Time", "Longitude"),
#'      ylab = c("Normalized Temperature (\u00B0C)", "Latitude"))
#'
#' plotly_funts(Montana[1:100],
#'              xlab = c("Time", "Longitude"),
#'              ylab = c("Normalized Temperature (\u00B0C)", "Latitude"),
#'              zlab = c("", "NDVI"),
#'              main = c("Temperature Curves", "NDVI Images"),
#'              color_palette = "RdYlGn",
#'              xticklabels = list(c("00:00","06:00","12:00","18:00","24:00"), c("113.40\u00B0 W", "113.30\u00B0 W")),
#'              xticklocs =list(c(1,6,12,18,24),c(1,33)),
#'              yticklabels = list(NA,c("48.70\u00B0 N", "48.77\u00B0 N")),
#'              yticklocs =list(NA,c(1,33)))
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
#'
#' plotly_funts(U$Lsingf[[1]])
#' plot(U$Lsingf[[2]])
#'
#' # MFSSA Reconstruction step:
#' Q <- freconstruct(U = U, groups = list(1, 2, 3))
#' plotly_funts(Q[[1]],
#'              xlab = c("Time", "Longitude"),
#'              ylab = c("Normalized Temperature (\u00B0C)", "Latitude"),
#'              zlab = c("", "NDVI"),
#'              main = c("Temperature Curves Mean", "NDVI Images Mean"), color_palette = "RdYlGn",
#'              xticklabels = list(c("00:00","06:00","12:00","18:00","24:00"),
#'                                 c("113.40\u00B0 W", "113.30\u00B0 W")),xticklocs =
#'                list(c(1,6,12,18,24),c(1,33)),
#'              yticklabels = list(NA,c("48.70\u00B0 N", "48.77\u00B0 N")),yticklocs =
#'                list(NA,c(1,33))) # mean
#'
#' plotly_funts(Q[[2]],
#'              xlab = c("Time", "Longitude"),
#'              ylab = c("Normalized Temperature (\u00B0C)", "Latitude"),
#'              zlab = c("", "NDVI"),
#'              main = c("Temperature Curves Periodic", "NDVI Images Periodic"), color_palette = "RdYlGn",
#'              xticklabels = list(c("00:00","06:00","12:00","18:00","24:00"),
#'                                 c("113.40\u00B0 W", "113.30\u00B0 W")),xticklocs =
#'                list(c(1,6,12,18,24),c(1,33)),
#'              yticklabels = list(NA,c("48.70\u00B0 N", "48.77\u00B0 N")),yticklocs =
#'                list(NA,c(1,33))) # periodic
#'
#' plot(Q[[3]], obs = 3,
#'      xlab = c("Time", "Longitude"),
#'      ylab = c("Normalized Temperature (\u00B0C)", "Latitude"),
#'      main = c("Temperature Curves Trend", "NDVI Images Trend,")) # trend
#' }
#' @useDynLib Rfssa
#' @export
fssa <- function(Y, L = NA, ntriples = 20, type = "ufssa") {
  N <- Y$N
  if (is.na(L)) L <- floor(N / 2L)
  if (ntriples > L){
    ntriples = L
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
    if(Y$dimSupp[[1]] == 1) {
      X_mat <- out[[i]]
    } else {
      n1 <- length(unique(Y$argval[[1]][,1]))
      n2 <- length(unique(Y$argval[[1]][,2]))
      X_mat <- array(out[[i]], dim = c(n1,n2,L))
    }
    Lsingf[[i]] <- funts(X = X_mat, basisobj = Y$basis[[1]])
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
    shifter[2L, j + 1L] <- shifter[2L, j] + L*ncol(Y$B_mat[[j]])
    Y_d[1,j] = ncol(Y$B_mat[[j]])
  }
  d <- cbind(0,Y_d*L)
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
      if(Y$dimSupp[[j]] == 1) {
        X_mat[[j]] <- my_pcs[[j]]
      } else {
        n1 <- length(unique(Y$argval[[j]][,1]))
        n2 <- length(unique(Y$argval[[j]][,2]))
        X_mat[[j]] <- array(my_pcs[[j]], dim = c(n1,n2,L))
      }
    }
    Lsingf[[i]] <- funts(X = X_mat, basisobj = Y$basis)
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
