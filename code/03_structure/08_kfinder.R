###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: run kfinder for all runs
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
#read file with parameters
#1. create parameters file
# vector with folder names for each run
runs <- dir("data/final", pattern = "run_", full.names = T)

runs %>%
  sapply(function(x){
    #get k values
    ks <-
      dir(x) %>%
        stringr::str_extract(pattern = "K[0-9]+") %>%
          stringr::str_extract(pattern = "[0-9]+") %>%
          as.numeric %>% range(na.rm = T)
    #get number of reps
    reps <-
      dir(x) %>%
        stringr::str_extract(pattern = "rep[0-9]+") %>%
          stringr::str_extract(pattern = "[0-9]+") %>%
          as.numeric %>% max(na.rm = T)
    #create parameter file
    param <- c(sprintf('"%s"', file.path(getwd(), x)),
               sprintf('"%s"', "str*"),
               ks[1], #from K
               ks[2], #to K
               readLines(file.path(x, "K2_rep1.stlog")) %>%
                 grep(pattern = " individuals$", value = T) %>%
                 stringr::str_extract("[0-9]+") %>%
                 as.numeric(),
               reps) #number of replicates per K
    fp1 <- file.path(x, "kfinder.par") #file path of parameter file
    writeLines(con = fp1, text = param) #write parameter file
#2. run kfinder within each run folder
system(paste("kfinder.mac", fp1))
  })
