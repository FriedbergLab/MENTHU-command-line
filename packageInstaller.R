# create local user library path (not present by default)
dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)

#packages <- c("stringr", "stringi", "Biostrings", "rentrez", "rlist", "plyr", "Rcpp", "curl", "httr", "jsonlite", "xml2")
packages <- c("stringr", "stringi", "BiocManager", "rentrez", "rlist", "plyr", "Rcpp", "curl", "httr", "jsonlite", "xml2")


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], lib = Sys.getenv("R_LIBS_USER"), repos = "http://cran.us.r-project.org")
}
BiocManager::install("Biostrings", lib = Sys.getenv("R_LIBS_USER"))

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
