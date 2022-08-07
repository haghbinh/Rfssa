#--------------------------------------------------------------
#' Functional Time Series Evaluation
#'
#' This is a method that can be used to evaluate a functional time series (\code{\link{fts}}) object at
#' new, specified grid points.
#'
#' @return A list of length p of matrices or three-dimensional arrays where each list entry type depends on whether the variable is observed over a one or two-dimensional domain.
#' @param Y An object of class \code{\link{fts}} to be evaluated.
#' @param grid A list of length p where each entry is a numeric, matrix, list, or \code{NULL} specifying the grid to evaluate the variable over. If the entry is \code{NULL}, then the variable is evaluated over the original grid.
#' @param ud_basis A list of length p where each entry is a matrix specifying the user-defined basis to be used for the evaluation of a variable. Note that the specified basis should match with the basis used to estimate the \code{\link{fts}} object in type and dimension.
#'
#' @examples
#' \dontrun{
#' ## Evaluate a univariate FTS at less grid points
#' require(Rfssa)
#' load_github_data("https://github.com/haghbinh/Rfssa/blob/master/data/Callcenter.RData")
#' ## Define functional objects
#' D <- matrix(sqrt(Callcenter$calls), nrow = 240)
#' N <- ncol(D)
#' time <- substr(seq(ISOdate(1999, 1, 1), ISOdate(1999, 12, 31), by = "day"),1,10)
#' K <- nrow(D)
#' u <- seq(0, K, length.out = K)
#' d <- 22 # Optimal Number of basis elements
#' ## Define functional time series
#' Y <- Rfssa::fts(list(D), list(list(d, "bspline")), list(u),time)
#' plot(Y, mains = c("Call Center Data Line Plot"),
#' xlabels = "Time (6 minutes aggregated)",
#' ylabels = "Sqrt of Call Numbers",type="line",
#' xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),xticklocs =
#'  list(c(1,60,120,180,240)))
#' u <- seq(0,K,length.out = 48)
#' D = eval.fts(Y=Y,grid = list(u))
#' Y <- Rfssa::fts(D, list(list(d, "bspline")), list(u))
#' plot(Y, mains = c("Call Center Data Line Plot"),
#' xlabels = "Time (6 minutes aggregated)",
#' ylabels = "Sqrt of Call Numbers",type="line",
#' xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),xticklocs =
#'  list(c(1,12,24,36,48)))
#'
#' ## Evaluate a multivariate FTS at more grid points
#'
#' require(Rfssa)
#' load_github_data("https://github.com/haghbinh/Rfssa/blob/master/data/Montana.RData")
#' Temp <- Montana$Temp
#' NDVI <- Montana$NDVI
#' d_temp <- 11
#' d_NDVI <- 13
#' ## Define functional time series
#' Y <- Rfssa::fts(
#'   list(Temp / sd(Temp), NDVI), list(
#'     list(d_temp, "bspline"),
#'     list(d_NDVI, d_NDVI, "bspline", "bspline")
#'   ),
#'   list(c(0, 23), list(c(1, 33), c(1, 33))),
#' time=colnames(Temp))
#' # Plot the first 100 observations
#' plot(Y[1:100],
#'      xlabels = c("Time", "Longitude"),
#'      ylabels = c("Standardized Temperature (\u00B0C)", "Latitude"),
#'      zlabels = c("", "NDVI"),
#'      mains = c("Temperature Curves", "NDVI Images"), color_palette = "RdYlGn",
#'      xticklabels = list(c("00:00","06:00","12:00","18:00","24:00"),
#'      c("113.40\u00B0 W", "113.30\u00B0 W")),xticklocs =
#'        list(c(1,6,12,18,24),c(1,33)),
#'        yticklabels = list(NA,c("48.70\u00B0 N", "48.77\u00B0 N")),yticklocs =
#'        list(NA,c(1,33))
#' )
#'
#' grid_temp = seq(0,23,length.out = 100)
#' u_NDVI = seq(1,33,length.out = 100)
#' v_NDVI = seq(1,33,length.out = 100)
#' grid_NDVI = list(u_NDVI,v_NDVI)
#' D = eval.fts(Y = Y, grid = list(grid_temp, grid_NDVI))
#'
#' Y <- Rfssa::fts(
#'   D, list(
#'     list(d_temp, "bspline"),
#'     list(d_NDVI,d_NDVI,"bspline","bspline")
#'   ),
#'   list(grid_temp, grid_NDVI)
#' )
#'
#' plot(Y[1:100],
#'      xlabels = c("Time", "Longitude"),
#'      ylabels = c("Standardized Temperature (\u00B0C)", "Latitude"),
#'      zlabels = c("", "NDVI"),
#'      mains = c("Temperature Curves", "NDVI Images"), color_palette = "RdYlGn",
#'      xticklabels = list(c("00:00","06:00","12:00","18:00","24:00"),
#'      c("113.40\u00B0 W", "113.30\u00B0 W")),xticklocs =
#'        list(c(1,24,48,72,100),c(1,33)),
#'        yticklabels = list(NA,c("48.70\u00B0 N", "48.77\u00B0 N")),yticklocs =
#'        list(NA,c(1,33))
#' )
#'
#' }
#'
#' @importFrom fda create.bspline.basis create.fourier.basis eval.basis
#' @export
eval.fts <- function(Y, grid, ud_basis = NULL){

  p = length(Y@C)
  Y_eval = list()
  basis = Y@basis_type
  if(is.null(ud_basis)) ud_basis = rep(NA,p);
  for(j in 1:p){
    C_j = Y@C[[j]]
    B_j = Y@B[[j]]

    if(basis[[j]][[1]]=="User-defined basis." && is.na(ud_basis[[j]][[1]])){

      stop(paste0("The basis for variable ", as.character(j)," needs to be provided by the user."))

    }else if(basis[[j]][[1]]=="User-defined basis." && is.na(ud_basis[[j]][[1]])==FALSE){

      Y_eval[[j]] = ud_basis[[j]]%*%C_j


    }else{

      if(ncol(Y@grid[[j]])==2){

        if(basis[[j]][[3]]=="bspline" && basis[[j]][[4]]=="bspline"){

          basis_new_1 = eval.basis(evalarg = grid[[j]][[1]],basisobj = create.bspline.basis(rangeval = c(min(Y@grid[[j]][,1]),max(Y@grid[[j]][,1])),nbasis = basis[[j]][[1]]))
          basis_new_2 = eval.basis(evalarg = grid[[j]][[2]],basisobj = create.bspline.basis(rangeval = c(min(Y@grid[[j]][,2]),max(Y@grid[[j]][,2])),nbasis = basis[[j]][[2]]))
          basis_new = kronecker(basis_new_1, basis_new_2)

        }else if(basis[[j]][[3]]=="fourier" && basis[[j]][[4]]=="bspline"){

          basis_new_1 = eval.basis(evalarg = grid[[j]][[1]],basisobj = create.fourier.basis(rangeval = c(min(Y@grid[[j]][,1]),max(Y@grid[[j]][,1])),nbasis = basis[[j]][[1]]))
          basis_new_2 = eval.basis(evalarg = grid[[j]][[2]],basisobj = create.bspline.basis(rangeval = c(min(Y@grid[[j]][,2]),max(Y@grid[[j]][,2])),nbasis = basis[[j]][[2]]))
          basis_new = kronecker(basis_new_1, basis_new_2)

        }else if(basis[[j]][[3]]=="bspline" && basis[[j]][[4]]=="fourier"){

          basis_new_1 = eval.basis(evalarg = grid[[j]][[1]],basisobj = create.bspline.basis(rangeval = c(min(Y@grid[[j]][,1]),max(Y@grid[[j]][,1])),nbasis = basis[[j]][[1]]))
          basis_new_2 = eval.basis(evalarg = grid[[j]][[2]],basisobj = create.fourier.basis(rangeval = c(min(Y@grid[[j]][,2]),max(Y@grid[[j]][,2])),nbasis = basis[[j]][[2]]))
          basis_new = kronecker(basis_new_1, basis_new_2)

        }else if(basis[[j]][[3]]=="fourier" && basis[[j]][[4]]=="fourier"){

          basis_new_1 = eval.basis(evalarg = grid[[j]][[1]],basisobj = create.fourier.basis(rangeval = c(min(Y@grid[[j]][,1]),max(Y@grid[[j]][,1])),nbasis = basis[[j]][[1]]))
          basis_new_2 = eval.basis(evalarg = grid[[j]][[2]],basisobj = create.fourier.basis(rangeval = c(min(Y@grid[[j]][,2]),max(Y@grid[[j]][,2])),nbasis = basis[[j]][[2]]))
          basis_new = kronecker(basis_new_1, basis_new_2)

          }

        Y_eval_j = array(data = basis_new%*%C_j, dim = c(length(grid[[j]][[1]]),length(grid[[j]][[2]]),ncol(C_j)))
        rotate <- function(x) t(apply(x, 2, rev))
        for(i in 1:dim(Y_eval_j)[3]) Y_eval_j[,,i]=rotate(Y_eval_j[,,i])[,c(ncol(Y_eval_j[,,i]):1)];
        Y_eval[[j]] = Y_eval_j
      }else if(ncol(Y@grid[[j]])==1){

        if(basis[[j]][[2]]=="bspline"){

          basis_new = eval.basis(evalarg = grid[[j]],basisobj = create.bspline.basis(rangeval = c(min(Y@grid[[j]]),max(Y@grid[[j]])),nbasis = basis[[j]][[1]]))

        }else if(basis[[j]][[2]]=="fourier"){

          basis_new = eval.basis(evalarg = grid[[j]],basisobj = create.fourier.basis(rangeval = c(min(Y@grid[[j]]),max(Y@grid[[j]])),nbasis = basis[[j]][[1]]))

        }
        Y_eval[[j]] = basis_new%*%C_j

      }else{
        Y_eval[[j]] = B_j%*%C_j
      }
    }

  }
  return(Y_eval)
}
