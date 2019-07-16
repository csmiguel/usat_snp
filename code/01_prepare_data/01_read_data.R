###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: read dartseq and usat genotypes and metadata and save them as R objects
#DESCRIPTION:
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
#  REQUIRED FILES:
#   Description: read data/raw/readme
#   Inpath: data/raw/
#  OUTPUT:
#    Description: raw genotypes plus metadata for dartseq samples
#    Outpath: data/intermediate/raw_genotypes.rds & metadata_dartseq.rds
#  DEPENDENCIES: functions/read_metadata.r & read_usats.r
###.............................................................................

library(dartR)
library(dplyr)
library(magrittr)

#1. DArTseq
#path to singlerow dartseq SNPs genotypes
path_2_genotypes <- c(
  "data/raw/Report_DHyl17-2984_SNP_mapping_2.csv",
  "data/raw/Report_DPelo19-4013_SNP_mapping_2.csv")

#read singlerow dartseq file into genlight object.
gen_dart <- lapply(path_2_genotypes, dartR::gl.read.dart)
names(gen_dart) <- c("hyla", "pelo")

#remove replica from Pelobates:
    #for Pelobates Ana included a blind
    #replica named 108541 from sample 106854.
gen_dart$pelo <- gl.drop.ind(gen_dart$pelo, "108541")
#read metadata
source("code/functions/read_metadata.r")
gen_dart %<>% metadata_Dartseq_Hyla() %>% metadata_Dartseq_Pelobates()
#2. microsatellite
source("code/functions/read_usats.r")
  #genlight objects with metadata
hyla_usat <- read_Hyla_molleri() #outputs metadata_hyla to global env
pelo_usat <- read_Pelobates_cultripes()#outputs metadata_pelo to global env

#create objects for exporting
gen <- list(dart_hyla = gen_dart[[1]], dart_pelo = gen_dart[[2]],
   usat_hyla = hyla_usat, usat_pelo = pelo_usat)
metadata_dartseq <- list(metadata_darthyla = metadata_darthyla,
  metadata_dartpelo = metadata_dartpelo)
#note: usat metadata is already stored inside genind objects.

saveRDS(gen, file = "data/intermediate/raw_genotypes.rds")
saveRDS(metadata_dartseq, file = "data/intermediate/metadata_dartseq.rds")
