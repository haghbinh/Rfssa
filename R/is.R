#' Check if an object is of class 'funts'
#'
#' @param obj The object to check.
#' @return TRUE if the object is of class 'funts', FALSE otherwise.
#'
#' @examples
#' \dontrun{
#' data("Callcenter")
#' is.funts(Callcenter)
#' }
#' @export
is.funts <- function(obj) {
  inherits(obj, "funts")
}
