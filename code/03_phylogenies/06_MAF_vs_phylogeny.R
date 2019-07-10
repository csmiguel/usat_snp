###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: analysis of bootstrap support
#DESCRIPTION:
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(adegenet)
library(dplyr)
library(ape)
library(dartR)

input_path <- "data/intermediate/gen_consolidated_filtered.rds"
gen <- readRDS(file = input_path)

source("code/parameters/boot.r")

#convert from genind to genlight
pelo <- gen[["dart_pelo"]] %>%
  dartR::gi2gl

#filter pelobates genotypes accross MAF thresholds
pelo_maf <-
  seq_along(maf) %>%
    lapply(function(x){
      dartR::gl.filter.maf(pelo, threshold = maf[x], v = 5)
      })

names(pelo_maf) <- as.character(maf)
