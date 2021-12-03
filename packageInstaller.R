# create local user library path (if not present by default)
if (!dir.exists(Sys.getenv('R_LIBS_USER'))){
  dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
} else {
  paste("Installing packages in local directory", Sys.getenv("R_LIBS_USER"))
}

# List required packages
#packages <- c("stringr", "stringi", "Biostrings", "rentrez", "rlist", "plyr", "Rcpp", "curl", "httr", "jsonlite", "xml2")
packages <- c("stringr", "stringi", "BiocManager", "rentrez", "rlist", "plyr", "Rcpp", "curl", "httr", "jsonlite", "xml2")


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], lib = Sys.getenv("R_LIBS_USER"))
}
BiocManager::install("Biostrings", lib = Sys.getenv("R_LIBS_USER"))

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
