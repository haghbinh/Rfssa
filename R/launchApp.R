#' Launch the Shiny Application Demonstration
#'
#' This function launches an app that can be used to help an individual better understand univariate or multivariate functional
#' singular spectrum analysis (\code{\link{fssa}}). The app allows the user to run univariate or multivariate functional singular spectrum analysis (depending on
#' the entered type of parameter)
#' on a variety of data types including simulated and real data available through the server. The app also has functionality
#' that allows the user to upload their own data. The app allows the user to compare different methods simultaneously such as
#' multivariate singular spectrum analysis versus univariate functional singular spectrum analysis. It also allows the user to choose
#' the number and types of basis elements used to estimate functional time series (\code{\link{fts}}) objects. The app supports \code{\link{fts}} plots and \code{\link{fssa}}
#' plots.
#'
#' @export launchApp
#'
#' @return A shiny application object.
#' @param type Type of FSSA with options of \code{type = "ufssa"} or \code{type = "mfssa"}.
#' @importFrom fda fd inprod eval.fd smooth.basis is.fd
#' @examples
#' \dontrun{
#'
#' launchApp()
#' }
#'
#' @importFrom plotly renderPlotly plotlyOutput
#' @importFrom fda pca.fd eval.fd
#' @importFrom graphics abline axis par plot points polygon
#' @importFrom stats fft integrate rnorm sd ts.plot density
#' @importFrom utils head read.table
#' @importFrom markdown mark
#' @importFrom Rssa ssa reconstruct
#' @import shiny
#'

# wrapper for shiny::shinyApp()
launchApp <- function(type = "ufssa") {
  options(shiny.sanitize.errors = FALSE)
  if (type == "ufssa") {
    shinyApp(ui = ui.fssa, server = server.fssa)
  } else if (type == "mfssa") {
    shinyApp(ui = ui.mfssa, server = server.mfssa)
  } else {
    stop("The type for the shiny is not valid.")
  }
}
