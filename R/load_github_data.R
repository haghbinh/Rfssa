#' Load Data from GitHub Repositories
#'
#' This function was found in
#' https://stackoverflow.com/questions/24846120/importing-data-into-r-rdata-from-github and can be
#' used to load .RData files from GitHub repositories. This function can be used to load the Callcenter,
#' Jambi, and Montana datasets from the Rfssa package hosted by GitHub at https://github.com/haghbinh/Rfssa.
#' It was found that hosting such datasets on GitHub and not including them in the Rfssa R package saved
#' a significant amount of space.
#'
#' @export load_github_data
#'
#' @return A dataset specified by a GitHub url.
#' @param github_data_url The GitHub url of the dataset.
#' @examples
#' \dontrun{
#' # Loading different datasets from the Rfssa repository hosted by GitHub.
#' call <- load_github_data("https://github.com/haghbinh/Rfssa/blob/master/data/Callcenter.RData")
#' jambi <- load_github_data("https://github.com/haghbinh/Rfssa/blob/master/data/Jambi.RData")
#' montana <- load_github_data("https://github.com/haghbinh/Rfssa/blob/master/data/Montana.RData")
#' }
#'
#' @importFrom httr GET stop_for_status content timeout
#'

load_github_data <- function(github_data_url) {
  url_len <- nchar(github_data_url)
  if (substr(github_data_url, start = (url_len - 8), stop = url_len) != "?raw=true") github_data_url <- paste0(github_data_url, "?raw=true")
  # based very closely on code for devtools::source_url
  temp_file <- tempfile()
  on.exit(unlink(temp_file))
  request <- httr::GET(github_data_url, httr::timeout(30))
  httr::stop_for_status(request)
  writeBin(httr::content(request, type = "raw"), temp_file)
  load(temp_file, envir = .GlobalEnv)
}
