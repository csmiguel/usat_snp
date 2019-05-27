###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: export data into structure format
#DESCRIPTION: export data into structure format and create params files
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
#  REQUIRED FILES:
#   Description:
#   Inpath:
#  OUTPUT:
#    Description: raw genotypes plus metadata for dartseq samples
#    Outpath:
#  DEPENDENCIES:
###.............................................................................
library(dartR)
library(dplyr)

input_path <- "data/intermediate/gen_consolidated_filtered.rds"
gen <- readRDS(file = input_path)
source("code/functions/edit_params.r")

#1. Create STRUCTURE input files from dart and usat genotypes
create_str_input(mode = "normal")

#2. Create STRUCTURE input files from SNP downsampling genotypes
  #genotypes for pelobates
pelo <- gen$dart_pelo
  #load subsampling sizes
source("code/parameters/boot.r")
s <- s[!(s %in% c(15000, 25000))]#edit vector with subsampling

create_str_input_subsets(mode = "normal")
