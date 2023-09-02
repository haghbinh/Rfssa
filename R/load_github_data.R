#' Load Jambi Data from GitHub Repository
#'
#' This function is used to load the Jambi dataset from the Rfssa package hosted on GitHub at
#' https://github.com/haghbinh/Rfssa. Hosting datasets on GitHub instead of including them
#' in the Rfssa R package saves a significant amount of space.
#'
#' @examples
#' \dontrun{
#' loadJambiData()
#' }
#'
#'
#' @export
loadJambiData <- function() {
  # Check if the data file exists locally
  if (file.exists("Jambi.RData")) {
    load("Jambi.RData", envir = .GlobalEnv)
  } else {
    # Download the data from GitHub
    url <- "https://github.com/haghbinh/Rfssa/raw/master/data/Jambi.RData"
    download.file(url, "Jambi.RData")
    load("Jambi.RData", envir = .GlobalEnv)
    save(Jambi, file = "Jambi.RData")
  }
}


#' Load Montana Data from GitHub Repository
#'
#' This function is used to load the Montana dataset from the Rfssa package hosted on GitHub at
#' https://github.com/haghbinh/Rfssa. Hosting datasets on GitHub instead of including them
#' in the Rfssa R package saves a significant amount of space.
#'
#' @examples
#' \dontrun{
#' loadMontanaData()
#' }
#'
#' @export
loadMontanaData <- function() {
  # Check if the data file exists locally
  if (file.exists("Montana.RData")) {
    load("Montana.RData", envir = .GlobalEnv)
  } else {
    # Download the data from GitHub
    url <- "https://github.com/haghbinh/Rfssa/raw/master/data/Montana.RData"
    download.file(url, "Montana.RData")
    load("Montana.RData", envir = .GlobalEnv)
    save(Montana, file = "Montana.RData")
  }
}
