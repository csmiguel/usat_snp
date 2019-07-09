###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: collect clumpak output for major cluster
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(dplyr)
#Collect Clumpak output:
datasets <-
  dir("data/intermediate/clumpak", pattern = "run_", full.names = T) %>%
    {.[-grep(., pattern = ".zip|_file")]}

#create list with Major Cluster ancestries
qclumpak <-
  seq_along(datasets) %>%
  lapply(function(x){
    #kruns <- dir(path = datasets[x], pattern = "K=[0-9]")
    k <- list()
    for (i in 1:8){
      k[[i]] <-
        paste0(datasets[x], "/K=", i,
          "/MajorCluster/CLUMPP.files/ClumppIndFile.output") %>%
        readLines() %>%
        gsub(pattern = "^.*: ", replacement = "") %>%
        textConnection %>%
        read.table
      names(k)[i] <- paste0("K", i)
    }
    k
  })

names(qclumpak) <- datasets %>% gsub(pattern = "^.*\\/run_", replacement = "")

saveRDS(qclumpak, "data/intermediate/clumpak_major.rds")
