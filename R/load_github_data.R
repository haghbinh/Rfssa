#' Load Callcenter Data from GitHub Repository
#'
#' This function retrieves the Callcenter dataset from the Rfssa_dataset repository on GitHub
#' (https://github.com/haghbinh/dataset/Rfssa_dataset).
#' The Callcenter dataset represents a small call center for an anonymous bank.
#'  It provides precise call timing data from January 1 to December 31, 1999.
#'  The data is aggregated into 6-minute intervals on each day.
#'  The returned object is a raw dataset in dataframe format;
#'  it is not a 'funts' class object.
#'  This raw data can then be further processed and converted into a 'funts' object named 'Callcenter'.
#'  See \code{\link{funts}} for more details on
#'  working with functional time series of class 'funts'.
#'
#' @format a dataframe with 87,600 rows and 5 variables:
#' \describe{
#'   \item{calls}{number of calls in a 6-minute aggregated interval.}
#'   \item{u}{numeric vector indicating the aggregated interval.}
#'   \item{Date}{date and time of call count recording.}
#'   \item{Day}{weekday associated with Date.}
#'   \item{Month}{month associated with Date.}
#' }

#'
#' @references
#' \enumerate{
#' \item
#' Brown, L., Gans, N., Mandelbaum, A., Sakov, A., Shen, H., Zeltyn, S., & Zhao, L. (2005).
#' Statistical analysis of a telephone call center: A queueing-science perspective.
#' \emph{Journal of the American Statistical Association}, \strong{100}(469), 36-50.
#' \item
#' Shen, H., & Huang, J. Z. (2005).
#' Analysis of call center arrival data using singular value decomposition.
#' \emph{Applied Stochastic Models in Business and Industry}, \strong{21}(3), 251-263.
#' \item
#' Huang, J. Z., Shen, H., & Buja, A. (2008).
#' Functional principal components analysis via penalized rank one approximation.
#' \emph{Electronic Journal of Statistics}, \strong{2}, 678-695.
#' \item
#' Maadooliat, M., Huang, J. Z., & Hu, J. (2015).
#' Integrating data transformation in principal components analysis.
#' \emph{Journal of Computational and Graphical Statistics}, \strong{24}(1), 84-103.
#' }
#'
#' @examples
#' # Load Callcenter data
#' Call_data <- loadCallcenterData()
#' D <- matrix(sqrt(Call_data), nrow = 240)
#'
#' # Define basis functions
#' bs1 <- create.bspline.basis(c(0, 23), 22)
#'
# Convert to `funts` object
#' Y <- funts(X = D, basisobj = bs1)
#'
#' @seealso \code{\link{funts}}
#'
#' @export
loadCallcenterData <- function() {
  fil <- tempfile("Callcenter", fileext = ".rds")
  # Download the data from GitHub
  url <- "https://github.com/haghbinh/dataset/raw/main/Rfssa_dataset/Callcenter.rds"
  download.file(url, destfile = fil)
  Callcenter <- readRDS(fil)
  return(Callcenter)
}
# =======================================================================
#'
#' Load Jambi MODIS Data from GitHub Repository
#'
#' This function retrieves the `Jambi` dataset from a GitHub repository hosted at
#' https://github.com/haghbinh/dataset/Rfssa_dataset. The Jambi dataset contains normalized difference vegetation index (NDVI) and enhanced vegetation
#' index (EVI) image data from NASA’s MODerate-resolution Imaging Spectroradiometer (MODIS) with
#' global coverage at a 250 m^2 resolution. The dataset covers the Jambi Province, Indonesia,
#' known for various forested land uses, including natural forests and plantations.
#' Monitoring land cover changes is crucial, especially in the context of forest exploitation and
#' conservation efforts. Seasonal variations significantly impact long-term land cover changes.
#' Data collection began on February 18, 2000, and continued until July 28, 2019, with data recorded
#' every 16 days. This dataset is valuable for studying vegetative land cover changes in the region.
#' The returned object is a raw dataset in `list` format;
#'
#' @format A list containing two arrays, each with dimensions 33 by 33 by 448. One array represents NDVI
#' image data, and the other represents EVI image data. The list also contains a date vector of length 448,
#' specifying the capture date for each 33 by 33 image.
#'
#' @references
#' \enumerate{
#'   \item Lambin, E., Geist, H., Lepers, E. (1999).
#'   Dynamics of Land-Use and Land-Cover Change in Tropical Regions.
#'   \emph{Annual Review of Environment and Resources}, 205-244.
#' }
#'
#' @source [MODIS Product Information](https://lpdaac.usgs.gov/products/mod13q1v006/)
#'
#' @examples
#' Jambi_raw <- loadJambiData()
#' str(Jambi_raw)
#'
#' @seealso  - The dataset object loaded by this function.
#'
#' @export
loadJambiData <- function() {
  fil <- tempfile("Jambi", fileext = ".rds")
  # Download the data from GitHub
  url <- "https://github.com/haghbinh/dataset/raw/main/Rfssa_dataset/Jambi.rds"
  download.file(url, destfile = fil)
  Jambi <- readRDS(fil)
  return(Jambi)
}
# =======================================================================

#' Load Montana Data from GitHub Repository
#'
#' This function retrieves the Montana dataset from a GitHub repository hosted at
#' https://github.com/haghbinh/dataset/Rfssa_dataset. Hosting datasets on GitHub
#' rather than including them in the Rfssa R package conserves storage space.
#' The Montana dataset contains intraday hourly temperature curves measured in degrees Celsius
#' and normalized difference vegetation index (NDVI) image data. Both types of data are recorded
#' near Saint Mary, Montana, USA. The NDVI images cover a region located between longitudes of
#' 113.30 degrees West and 113.56 degrees West and latitudes of 48.71 degrees North and 48.78 degrees North.
#' For each recorded intraday temperature curve, an NDVI image was captured on the same day every
#' 16 days, starting from January 1, 2008, and ending on September 30, 2013.
#' The dataset is valuable for environmental analysis, especially in the context of studying the impact
#' of temperature changes on vegetation. Combining both temperature and NDVI data can reveal more informative
#' patterns and insights.  The returned object is a raw dataset in `list` format;
#' This raw data can then be further processed and converted into a 'funts' object named 'Montana'.
#'  See \code{\link{funts}} for more details on  working with functional time series of class 'funts'.
#'
#'
#' @format A list containing two components:
#' \describe{
#'   \item{Temperature Data}{A 24 by 133 matrix of discrete samplings of intraday hourly temperature curves.}
#'   \item{NDVI Images}{An array with dimensions 33 by 33 by 133, where each 33 by 33 slice represents an NDVI image.}
#' }
#'
#' @references
#' \enumerate{
#'   \item Diamond, H. J., Karl, T., Palecki, M. A., Baker, C. B., Bell, J. E., Leeper, R. D.,
#'      Easterling, D. R., Lawrimore, J. H., Meyers, T. P., Helfert, M. R., Goodge, G.,
#'      and Thorne, P.W. (2013). U.S. climate reference network after one decade of operations:
#'      status and assessment. [Read More](https://www.ncdc.noaa.gov/crn/qcdatasets.html).
#'      Last accessed April 2020.
#'   \item Tuck, S. L., Phillips, H. R., Hintzen, R. E., Scharlemann, J. P., Purvis, A., and
#'      Hudson, L. N. (2014). MODISTools – downloading and processing MODIS
#'      remotely sensed data in R. Ecology and Evolution, 4(24):4658–4668.
#' }
#'
#' @examples
#' # Load Montana data
#' montana_data <- loadMontanaData()
#'
#' # Extract variables
#' Temp <- montana_data$Temp
#' NDVI <- montana_data$NDVI
#'
#' # Create a list for Montana data
#' Montana_Data <- list(Temp / sd(Temp), NDVI)
#'
#' # Define basis functions
#' bs1 <- create.bspline.basis(c(0, 23), 11)
#' bs2 <- create.bspline.basis(c(1, 33), 13)
#' bs2d <- list(bs2, bs2)
#' bsmv <- list(bs1, bs2d)
#'
#' # Convert to funts object
#' Y <- funts(X = Montana_Data, basisobj = bsmv,
#'            start = as.Date("2008-01-01"),
#'            end = as.Date("2013-09-30"),
#'            vnames = c("Normalized Temperature (\u00B0C)" , "NDVI"),
#'            dnames = list("Time", c("Latitude", "Longitude")),
#'            tname = "Date"
#' )
#'
#' @seealso \code{\link{funts}}
#'
#' @export
loadMontanaData <- function() {
  fil <- tempfile("Montana", fileext = ".rds")
  # Download the data from GitHub
  url <- "https://github.com/haghbinh/dataset/raw/main/Rfssa_dataset/Montana.rds"
  download.file(url, destfile = fil)
  Montana <- readRDS(fil)
  return(Montana)
}


# =======================================================================

#' Load Austin Temperature Data from GitHub Repository
#'
#' This function retrieves the Austin Temperature dataset from a GitHub repository hosted at
#' https://github.com/haghbinh/dataset/Rfssa_dataset. Hosting datasets on GitHub
#' rather than including them in the Rfssa R package conserves storage space.
#' The Austin Temperature dataset contains intraday hourly temperature curves
#' measured in degrees Celsius from January 1973 to July 2023, recorded once per month.
#' The returned object is a raw dataset in `list` format;
#'
#' @format Containing two matrix objects:
#' \describe{
#'   \item{Austin Temperature Data}{A 24 by 607 matrix of discrete samplings of intraday hourly temperature curves, one day per month.}
#'   \item{Utqiagvik Temperature Data}{A 24 by 607 matrix of discrete samplings of intraday hourly temperature curves, one day per month.}
#' }
#'
#' @references
#'   Trinka, J., Haghbin, H., Shang, H., Maadooliat, M. (2023). Functional Time Series Forecasting: Functional Singular Spectrum Analysis Approaches. Stat, e621.
#'
#' @examples
#' Austin_raw <- loadAustinData()
#' str(Austin_raw)
#'
#'
#' @export
loadAustinData <- function() {
  fil <- tempfile("Austin", fileext = ".rds")
  # Download the data from GitHub
  url <- "https://github.com/haghbinh/dataset/raw/main/Rfssa_dataset/Austin.rds"
  download.file(url, destfile = fil)
  Austin <- readRDS(fil)
  return(Austin)
}

# =======================================================================



#' Load Utqiagvik Temperature Data from GitHub Repository
#'
#' This function retrieves the Utqiagvik Temperature dataset from a GitHub repository hosted at
#' https://github.com/haghbinh/dataset/Rfssa_dataset. Hosting datasets on GitHub
#' rather than including them in the Rfssa R package conserves storage space.
#' The Utqiagvik Temperature dataset contains intraday hourly temperature curves
#' measured in degrees Celsius from January 1973 to July 2023, recorded once per month.
#' The returned object is a raw dataset in `list` format;
#' @format Containing two matrix objects:
#' \describe{
#'   \item{Austin Temperature Data}{A 24 by 607 matrix of discrete samplings of intraday hourly temperature curves, one day per month.}
#'   \item{Utqiagvik Temperature Data}{A 24 by 607 matrix of discrete samplings of intraday hourly temperature curves, one day per month.}
#' }
#'
#' @references
#'   Trinka, J., Haghbin, H., Shang, H., Maadooliat, M. (2023). Functional Time Series Forecasting: Functional Singular Spectrum Analysis Approaches. Stat, e621.
#'
#' @examples
#' Utqiagvik_raw <- loadUtqiagvikData()
#' str(Utqiagvik_raw)
#'
#' @export
loadUtqiagvikData <- function() {
  fil <- tempfile("Utqiagvik", fileext = ".rds")
  # Download the data from GitHub
  url <- "https://github.com/haghbinh/dataset/raw/main/Rfssa_dataset/Utqiagvik.rds"
  download.file(url, destfile = fil)
  Utqiagvik <- readRDS(fil)
  return(Utqiagvik)
}
