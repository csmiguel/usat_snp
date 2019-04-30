###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: produce phylogenetic trees for microsatellite and snp data
#DESCRIPTION:
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
library(adegenet)
library(assertthat)

input_path <- paste0("data/intermediate/filt_genotypes.rds")
gen <- readRDS(file = input_path)

source("code/functions/bootstrap_trees.r")

trees_gen <- list()
sapply(gen, function(x){
  type <- class(x)
  sp <- names(gen)
  bootstrap_nj(gen, 200, type)
}



  #calculate bootstrap support
  bootnj_Da <- ape::countBipartitions(true_tree, zz) #issue2
