<<<<<<< HEAD
#' Check if an object is of class 'funts'
#'
#' @param obj The object to check.
#' @return TRUE if the object is of class 'funts', FALSE otherwise.
#' @export
is.funts <- function(obj) {
  inherits(obj, "funts")
=======
#' Check if an object is of class 'funfts'
#'
#' @param obj The object to check.
#' @return TRUE if the object is of class 'funfts', FALSE otherwise.
#' @export
is.funfts <- function(obj) {
  inherits(obj, "funfts")
>>>>>>> 7551633c78ae9dba6b403d65fa72dcd5691d31de
}
