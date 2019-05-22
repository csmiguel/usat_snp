###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: populations within genotypes to as the usat and SNP datasets to match
#the most.
#DESCRIPTION:
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
#  REQUIRED FILES:
#   Description: genlight objects with raw genotypes for Hyla and Pelobates
#   Inpath: data/intermediate/raw_genotypes.rds
#  OUTPUT:
#    Description: raw genotypes in genlight format
#    Outpath: data/intermediate/filtered_genotypes.rds
###.............................................................................
library(dartR)
library(dplyr)

input_path <- "data/intermediate/filt_genotypes.rds"
gen <- readRDS(file = input_path)
input_path2 <- "data/intermediate/metadata_all.rds"
meta <- readRDS(input_path2)
source("code/functions/map_distribution.r")

#check metadata only contains samples in filtered genotypes:
  #assert_that:
assertthat::assert_that(
  seq_along(gen) %>% sapply(function(x){
    #all ids in genotypes are present in metadata
  all(indNames(gen[[x]]) %in% meta[[x]]$sample_id) &
    #all ids unique
  all(duplicated(meta[[x]]) == F) &
    #vectors with names have same length
  adegenet::nInd(gen[[x]]) == length(meta[[x]]$sample_id)
  }) %>% all())

#Plot samples in map so as to decide which populations to keep
plot_samples(meta = meta,
  name_plot = "data/intermediate/pre_sample_fil_distribution.pdf")

#based on plotting samples on a map the
