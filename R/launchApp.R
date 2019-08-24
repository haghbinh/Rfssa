#' launches the shinyAppDemo app
#'
#' @export launchApp
#'
#' @return shiny application object
#'
#' @example \dontrun {launchApp()}
#'
#' @import shiny
#'


# wrapper for shiny::shinyApp()
launchApp <- function(type="FSSA") {
  if (type=="FSSA")
    shinyApp(ui = ui.fssa, server = server.fssa)
  else
    shinyApp(ui = ui.mfssa, server = server.mfssa)
}
