###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: read kfinder results r
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

library(dplyr)
library(magrittr)
library(starmie)

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

# tidy kfinder results
kfinder_tidy <-
  kfinder_ls %>%
    plyr::ldply() %>%
    dplyr::as_tibble() %>%
    dplyr::rename(dataset = .id) %>%
    dplyr::mutate(marker = replace(dataset,
        stringr::str_detect(dataset, "dart"), "SNPs") %>%
        replace(stringr::str_detect(., "usat"), "microsatellites") %>%
        replace(stringr::str_detect(., "[0-9]"), "SNPs") %>% as.factor) %>%
    dplyr::mutate(species = replace(dataset,
        stringr::str_detect(dataset, "pelo"), "P. cultripes") %>%
        replace(stringr::str_detect(., "hyla"), "H. molleri") %>% as.factor) %>%
    dplyr::mutate(Evanno = suppressWarnings(as.numeric(as.character(Evanno))))

#2. best FILES
best_files <-
  kfinder %>%
  sapply(function(x){
    readLines(x) %>%
      grep(pattern = "Best output", value = T) %>%
      sub(pattern = "^.* /", replacement = "/", .)
  })

#3. Create Evanno plots and object to save
# read runs
runs <- dir("data/final", pattern = "run_", full.names = T)
# warning! very time consuming step:
h <-
  seq_along(runs) %>%
    lapply(function(x){
    #list files with runs within each run folder
    kreps <- list.files(runs[x], pattern = "_f$", full.names = TRUE)
    # read structure files into starmie object
    hh <- starmie::structList(lapply(kreps, starmie::loadStructure))
    # compute Evanno's metrics
    zz <- starmie::bestK(hh)
    # save plots
    ggplot2::ggsave(file.path(runs[x], "evanno.pdf"), device = "pdf")
    #return table with values
    zz
  })
names(h) <- sapply(runs, basename)

#tidy Evanno
evanno_tidy <-
  h %>%
  plyr::ldply() %>%
  dplyr::as_tibble() %>%
  dplyr::rename(dataset = .id) %>%
  dplyr::mutate(marker = replace(dataset,
            stringr::str_detect(dataset, "dart"), "SNPs") %>%
            replace(stringr::str_detect(., "usat"), "microsatellites") %>%
            replace(stringr::str_detect(., "[0-9]"), "SNPs") %>% as.factor) %>%
  dplyr::mutate(species = replace(dataset,
            stringr::str_detect(dataset, "pelo"), "P. cultripes") %>%
            replace(stringr::str_detect(., "hyla"), "H. molleri") %>% as.factor)

#save objects
saveRDS(kfinder_ls, "data/intermediate/kfinder_ls.rds")
saveRDS(kfinder_tidy, "data/intermediate/kfinder_tidy.rds")
saveRDS(best_files, "data/intermediate/best_files_kfinder.rds")
saveRDS(h, "data/intermediate/evanno.rds")
saveRDS(evanno_tidy, "data/intermediate/evanno_tidy.rds")
