#' Load packages and install them if necessary from BiocManager
#'
#' This function loads the provided list of packages, and installs them if necessary from BiocManager
#'
#' @param package1 Provide list of packages to load or install
#' @return List of packages that are loaded or installed
#' @import BiocManager
#' @import dplyr
#' @export

# Function to install&load/load required R packages
install_load <- function (package1, ...) {
  packages <- c(package1, ...)
  for(package in packages){
    if(package %in% rownames(installed.packages()))
      next
    else {
      install.packages(package)
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(package)
    }
  }
  lapply(packages, require, character.only = TRUE)
}

