#' Load Callcenter Data from GitHub Repository
#'
#' This function loads the Callcenter dataset from the Rfssa_dataset repository hosted on GitHub
#' (https://github.com/haghbinh/dataset/Rfssa_dataset). Hosting datasets on GitHub instead of
#' including them in the Rfssa R package saves a significant amount of space.
#' The Callcenter dataset represents a small call center for an anonymous bank, providing the exact time
#' of calls connected to the center from January 1 to December 31 in the year 1999. The data is aggregated
#' into time intervals, with each element of the data matrix representing call volume during a specific
#' 6-minute interval on a particular day.
#'
#' @format A dataframe with 87,600 rows and 5 variables:
#' \describe{
#'   \item{calls}{The number of calls in a 6-minute aggregated interval.}
#'   \item{u}{A numeric vector indicating the aggregated interval.}
#'   \item{Date}{Date and time when the call counts were recorded.}
#'   \item{Day}{Weekday associated with Date.}
#'   \item{Month}{Month associated with Date.}
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
#' \dontrun{
#' loadCallcenterData()
#' }
#'
#' @seealso \code{\link{Callcenter}}, \code{\link{funts}}
#'
#' @export
loadCallcenterData <- function() {
  # Check if the data file exists locally
  if (file.exists("Callcenter.RData")) {
    load("Callcenter.RData", envir = .GlobalEnv)
  } else {
    # Download the data from GitHub
    url <- "https://github.com/haghbinh/dataset/raw/main/Rfssa_dataset/Callcenter.RData"
    download.file(url, "Callcenter.RData")
    load("Callcenter.RData", envir = .GlobalEnv)
    save(Callcenter, file = "Callcenter.RData")
  }
}
# =======================================================================
#'
#' Load Jambi MODIS Data from GitHub Repository
#'
#' This function is used to load the Jambi dataset from a GitHub repository hosted at
#' https://github.com/haghbinh/dataset/Rfssa_dataset. By hosting datasets on GitHub
#' instead of including them in the Rfssa R package, it saves a significant amount of space.
#' The Jambi dataset contains normalized difference vegetation index (NDVI) and enhanced vegetation
#' index (EVI) image data from NASA’s MODerate-resolution Imaging Spectroradiometer (MODIS) with
#' global coverage at 250 m^2 resolution. The dataset covers the Jambi Province, Indonesia,
#' which is known for various forested land utilizations, including natural forests and plantations.
#' Monitoring such land cover changes is essential, especially in the context of forest exploitation and
#' conservation efforts. Seasonal variations play a significant role in long-term land cover changes.
#'
#' The data collection started on February 18, 2000, and continued until July 28, 2019, with data recorded
#' every 16 days. This dataset is valuable for studying vegetative land cover changes in the region.
#'
#' @format A list containing two arrays, each with dimensions 33 by 33 by 448. One array represents NDVI
#' image data, and the other represents EVI image data. The list also contains a date vector of length 448,
#' specifying the date of capture for each 33 by 33 image.
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
#' \dontrun{
#'   loadJambiData()
#' }
#'
#' @seealso \code{\link{Jambi}} - The dataset object loaded by this function.
#'
#' @export
loadJambiData <- function() {
  # Check if the data file exists locally
  if (file.exists("Jambi.RData")) {
    load("Jambi.RData", envir = .GlobalEnv)
  } else {
    # Download the data from GitHub
    url <- "https://github.com/haghbinh/dataset/raw/main/Rfssa_dataset/Jambi.RData"
    download.file(url, "Jambi.RData")
    load("Jambi.RData", envir = .GlobalEnv)
    save(Jambi, file = "Jambi.RData")
  }
}
# =======================================================================

#' Load Montana Data from GitHub Repository
#'
#' This function loads the Montana dataset from a GitHub repository hosted at
#' https://github.com/haghbinh/dataset/Rfssa_dataset. By hosting datasets on GitHub
#' instead of including them in the Rfssa R package, it saves a significant amount of space.
#'
#' The Montana dataset contains intraday hourly temperature curves measured in degrees Celsius
#' and normalized difference vegetation index (NDVI) image data. Both types of data are recorded
#' near Saint Mary, Montana, USA. The NDVI images cover a region located between longitudes of
#' 113.30 degrees West and 113.56 degrees West and latitudes of 48.71 degrees North and 48.78 degrees North.
#' For each recorded intraday temperature curve, an NDVI image was captured on the same day every
#' 16 days, starting from January 1, 2008, and ending on September 30, 2013.
#' The dataset is valuable for environmental analysis, especially in the context of studying the impact
#' of temperature changes on vegetation. Combining both temperature and NDVI data can reveal more informative
#' patterns and insights.
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
#' \dontrun{
#'   loadMontanaData()
#' }
#'
#' @seealso \code{\link{Montana}} - The dataset object loaded by this function.
#'
#' @export
loadMontanaData <- function() {
  # Check if the data file exists locally
  if (file.exists("Montana.RData")) {
    load("Montana.RData", envir = .GlobalEnv)
  } else {
    # Download the data from GitHub
    url <- "https://github.com/haghbinh/dataset/raw/main/Rfssa_dataset/Montana.RData"
    download.file(url, "Montana.RData")
    load("Montana.RData", envir = .GlobalEnv)
    save(Montana, file = "Montana.RData")
  }
}
