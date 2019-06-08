###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: calculate best K
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
#  REQUIRED FILES:
#   Description:
#   Inpath:
#  OUTPUT:
#    Description:
#    Outpath:
#  DEPENDENCIES:
###.............................................................................
library(dplyr)
#1. create parameters file
# vector with folder names for each run
kfinder <- dir("data/final", pattern = "str.K$", full.names = T, recursive = T)
#create list of data frames with K statistics
kfinder_ls <-
  kfinder %>%
    lapply(function(x){
      readLines(x) %>%
      grep(pattern = "^\\s+[0-9]", value = T) %>%
      textConnection() %>%
      read.table() %>%
      setNames(c("K", "Pritchard", "Evanno", "Parsimony"))
    })
names(kfinder_ls) <- kfinder %>%
                      stringr::str_extract("run.*/") %>%
                      gsub(pattern = "/", replacement = "")
#best FILES
best_files <-
  kfinder %>%
    sapply(function(x){
      readLines(x) %>%
        grep(pattern = "Best output", value = T) %>%
        sub(pattern = "^.* /", replacement = "/", .)
    })
