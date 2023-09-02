#' Load Data from GitHub Repositories
#'
#' These functions are used to load .RData files from GitHub repositories. They can be used to load the
#' Jambi and Montana datasets from the Rfssa package hosted on GitHub at https://github.com/haghbinh/Rfssa.
#' Hosting datasets on GitHub instead of including them in the Rfssa R package saves
#' a significant amount of space.
#'
#' @docType package
#'
#' @importFrom methods is new
#'
#' @examples
#' \dontrun{
#' loadJambiData()
#' loadMontanaData()
#' }

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

#' @export
loadMontanaData <- function() {
  # Check if the data file exists locally
  if (file.exists("Jambi.RData")) {
    load("Jambi.RData", envir = .GlobalEnv)
  } else {
    # Download the data from GitHub
    url <- "https://github.com/haghbinh/Rfssa/blob/master/data/Montana.RData"
    download.file(url, "Jambi.RData")
    load("Montana.RData", envir = .GlobalEnv)
    save(Montana, file = "Montana.RData")
  }
}
