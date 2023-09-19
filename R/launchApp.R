#' Launch Shiny Application for FSSA Demonstration
#'
#' This function launches a Shiny app to facilitate the understanding of univariate
#' or multivariate Functional Singular Spectrum Analysis (\code{\link{fssa}}). The app
#' enables users to perform univariate or multivariate FSSA on various data types,
#' including simulated and real data provided by the server. Users can also upload
#' their own data. The app supports simultaneous comparisons of different methods,
#' such as multivariate vs. univariate FSSA, and allows users to select the number
#' and types of basis elements used for estimating Functional Time Series
#' (\code{\link{funts}}) objects. It offers plotting capabilities for both
#'  \code{\link{funts}} and \code{\link{fssa}}.
#'
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
    stop("Invalid Shiny app type.")
  }
}
