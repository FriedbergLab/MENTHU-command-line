#packages <- c("stringr", "stringi", "Biostrings", "rentrez", "rlist", "plyr", "Rcpp", "curl", "httr", "jsonlite", "xml2")
#packages <- c("stringr", "stringi", "rentrez", "rlist", "plyr", "Rcpp", "curl", "httr", "jsonlite", "xml2")
packages <- c("stringr", "stringi", "BiocManager", "rentrez", "rlist", "plyr", "Rcpp", "curl", "httr", "jsonlite", "xml2")


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
BiocManager::install("Biostrings")

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("Biostrings")
