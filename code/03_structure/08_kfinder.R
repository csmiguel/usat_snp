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

#1. create parameters file
# vector with folder names for each run
runs <- dir("data/final", pattern = "run_", full.names = T)

runs %>%
  sapply(function(x){
    #create parameter file
    param <- c(sprintf('"%s"', file.path(getwd(), x)),
               sprintf('"%s"', "str*"),
               2, #from K
               8, #to K
               readLines(file.path(x, "K2_rep1.stlog")) %>%
                 grep(pattern = " individuals$", value = T) %>%
                 stringr::str_extract("[0-9]+") %>%
                 as.numeric(),
               5) #number of replicates per K
    fp1 <- file.path(x, "kfinder.par") #file path of parameter file
    writeLines(con = fp1, text = param) #write parameter file
#2. run kfinder within each run folder
system(paste("kfinder.mac", fp1))
  })
