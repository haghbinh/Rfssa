#' Jambi MODIS Data
#'
#' This data set contains the normalized difference vegetation index (NDVI) and enhanced vegetation index (EVI) image data from NASA’s MODerate-resolution Imaging Spectroradiometer (MODIS) with global coverage at 250 m^2.
#' The goal of the study is to collect raw image data of the Jambi Province, Indonesia.
#' Indonesia manages various forested land utilizations such as natural forest and plantations which, in the past, have been exploited throughout the country.
#' Greater criticisms on forest exploitation lead to a moratorium which needs to be
#' monitored frequently. Assessment of woody vegetation could be taken using field surveys or
#' remote sensing. It was found that season is probably the most intriguing factor in vegetative land cover,
#' especially in long-term land cover changes (Lambin, 1999).
#' The data was gathered starting in 2000-02-18 and ending in 2019-07-28 every 16 days.
#' @name Jambi
#' @format A list which contains two 33 by 33 by 448 arrays where one array is for NDVI image data and the other is for EVI image data. The list also contains a date vector of length 448 which specifies upon which date was each image 33 by 33 image taken.
#' \describe{
#'   \item{Days 1 - 448}{Pixel intensity with values between zero and one}
#' } @references
#' \enumerate{
#' \item
#'  Lambin, E., Geist, H., Lepers, E. (1999).
#'  Dynamics of Land-Use and Land-Cover Change in Tropical Regions
#'  \emph{Annual Review of Environment and Resources}, 205-244.
#'
#'
#' }
#' @source \url{https://lpdaac.usgs.gov/products/mod13q1v006/}
#' @seealso \code{\link{fssa}}
"Jambi"
