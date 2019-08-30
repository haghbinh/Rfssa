#' Jambi NDVI Density Function Data
#'
#' This data set is Normalized Difference Vegetation Index, NDVI, from NASA’s MODerate-resolution Imaging Spectroradiometer (MODIS) with global coverage at 250 m^2.
#' The study was located in Jambi Province, Indonesia.
#' Indonesia manages various forested land utilizations, for instance natural forest and plantations.
#' In the past, natural forest had been exploited throughout the country.
#' Greater criticisms on forest exploitation lead to a moratorium which needs to be
#' monitored frequently.Assessment of woody vegetation could be taken using field surveys or
#' remote sensing. Season is probably the most intriguing factor in vegetative land cover,
#' especially in long-term land cover changes (Lambin, 1999).
#' The data was gathered starting in 2000-01-01 and ending in 2019-05-12 every 16 days.
#' Therefore, the original dataset includes 441 NDVI images of size 33*33 pixels. Kernel desnity estimation was performed on this data to obtain 441 NDVI densities that give the distribution of pixel intensity in each image.
#' @name NDVI
#' @format A data matrix with 512 rows and 441 columns where each column is one of the days an image was taken:
#' \describe{
#'   \item{Days 1 - 441}{density of pixel intensity}
#' } @references
#' \enumerate{
#' \item
#'  Lambin, E., Geist, H., Lepers, E. (1999).
#'  Dynamics of Land-Use and Land-Cover Change in Tropical Regions
#'  \emph{Annual Review of Environment and Resources}, 205-244.
#'
#' }
"NDVI"
