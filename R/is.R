#' Check if an object is of class 'funts'
#'
#' @param obj The object to check.
#' @return TRUE if the object is of class 'funts', FALSE otherwise.
#' @export
is.funts <- function(obj) {
  inherits(obj, "funts")
}
