#' Number of calls for a bank.
#'
#' This dataset is a small call center for an anonymous bank (Brown et al.,
#' 2005). This dataset provides the exact time of the calls that were connected to
#' the call center between January 1 and December 31 in the year 1999.
#' The data are aggregated  into time intervals to obtain a data matrix. More
#' precisely, the i'th, row jth column of the data matrix contains the call volume
#' during the jth time interval on day i. This data set has been analyzed in several
#' prior studies; e.g. Brown et al. (2005), Huang et al. (2008), Shen and Huang
#' (2008) and Maadooliat et al. (2015). Here, the data are aggregated  into time
#' intervals 6 minutes, so a data matrix of dimension 240 by 365 is resulted.
<<<<<<< HEAD
#' @name Callcenter
#' @format A data frame with 87600 rows and 5 variables:
#' @name callcenter
#' @format A data frame with 240 rows and 366 variables:
=======
<<<<<<< HEAD
#' @name Callcenter
#' @format A data frame with 87600 rows and 5 variables:
=======
#' @name callcenter
#' @format A data frame with 240 rows and 366 variables:
>>>>>>> 39e05e79cd864f01dc2024c37a0786534d97e25a
>>>>>>> 1da5083a117ca26b5ef20f6607b189568dfe8481
#' \describe{
#'   \item{calls}{The number of calls in 6 minutes aggregated interval.}
#'   \item{u}{a numeric vector to show the aggregated interval}
#'   \item{Date}{Date time when the calls counts are recorded}
#'   \item{Day}{Weekday associated with Date}
#'   \item{Month}{Month associated with Date}
#' }
"Callcenter"
