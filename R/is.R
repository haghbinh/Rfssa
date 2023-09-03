#' Check if an object is of class 'funfts'
#'
#' @param obj The object to check.
#' @return TRUE if the object is of class 'funfts', FALSE otherwise.
#' @export
is.funfts <- function(obj) {
  inherits(obj, "funfts")
}
