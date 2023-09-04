#' Indexing into Functional Time Series
#'
#' An indexing method for functional time series (\code{\link{fts}}).
#' @return An object of class \code{\link{fts}}.
#'
#'
#' @param Y An object of class \code{\link{fts}}.
#' @param i The index value.
#' @examples
#' \dontrun{
#' require(Rfssa)
#' load_github_data("https://github.com/haghbinh/Rfssa/blob/master/data/Callcenter.RData")
#' D <- matrix(sqrt(Callcenter$calls), nrow = 240)
#' u <- seq(0, 1, length.out = 240) # Define domain of functional data
#' d <- 22 # number of basis elements
#' Y <- Rfssa::fts(list(D), list(list(d, "bspline")), list(u))
#' plot(Y)
#' plot(Y[10:15])
#' }
#'
#' @seealso \code{\link{fts}}
#' @useDynLib Rfssa
#' @export
"[.fts" <- function(Y, i = "index") {
  p <- length(Y@C)
  grid <- Y@grid
  basis <- Y@B
  eval_fts <- list()
  new_grid <- list()
  for (j in 1:p) {
    eval_fts[[j]] <- basis[[j]] %*% Y@C[[j]][, i]
    N <- ncol(eval_fts[[j]])

    if (ncol(Y@grid[[j]]) == 2) {
      x <- unique(Y@grid[[j]][, 1])
      y <- unique(Y@grid[[j]][, 2])
      new_fts_two_d <- array(data = NA, dim = c(length(x), length(y), N))
      for (n in 1:N) {
        count <- 1
        for (i_1 in 1:length(x)) {
          for (i_2 in 1:length(y)) {
            new_fts_two_d[i_1, i_2, n] <- eval_fts[[j]][count, n]
            count <- count + 1
          }
        }
      }

      eval_fts[[j]] <- new_fts_two_d
      new_grid[[j]] <- list(x, y)
    } else {
      new_grid[[j]] <- Y@grid[[j]]
    }
  }
  time = colnames(Y@C[[1]])[i]
  fts_out <- Rfssa::fts(eval_fts, basis, new_grid, time)
  fts_out@basis_type <- Y@basis_type
  out <- fts_out
  return(out)
}
