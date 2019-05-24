###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: consolidation of metadata
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
library(magrittr)
#filtered dartseq genotypes
input_path <- "data/intermediate/dart_filt.rds"
dart_gen <- readRDS(file = input_path)
#dartseq metadata
input_path2 <- "data/intermediate/metadata_dartseq.rds"
dart_meta <- readRDS(input_path2)
names(dart_meta) <- names(dart_gen)
#usat gen
input_path <- "data/intermediate/raw_genotypes.rds"
usat_gen <- readRDS(file = input_path) %>% {.[grep("usat", names(.))]}
#functions
source("code/functions/consolidate_metadata.r")

#replace @other in genlight by metadata
for (i in seq_along(dart_meta)){
other(dart_gen[[i]])$loc.metrics <- NULL
dart_gen[[i]]@other <- dart_meta[[i]]
}
rm(i, dart_meta)
#merge all genotypes objects:
gen <- c(dart_gen, usat_gen)
rm(dart_gen, usat_gen)

#CONSOLIDATE metadata: removes samples in metadata non present in
#indNames and reorders metadata$sample_id according to order in indNames
for (i in seq_along(gen)){
  gen[[i]] <- consolidate_metadata(gen[[i]])
}
