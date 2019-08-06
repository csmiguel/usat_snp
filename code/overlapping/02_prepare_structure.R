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

library(dartR)
library(dplyr)

input_path <- "data/intermediate/shared_gen_consolidated_filtered.rds"
gen <- readRDS(file = input_path)
names(gen) <- paste0("shared_", names(gen))
input_path2 <- "data/intermediate/lambda.rds"
lambda <- readRDS(input_path2)
source("code/functions/edit_params.r")
source("code/parameters/boot.r")

#1. Create STRUCTURE input files from dart and usat genotypes
create_str_input(mode = "normal")
