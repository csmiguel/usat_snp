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
library(dplyr)
library(ape)

input_path1 <- paste0("data/intermediate/bt_trees.rds")
input_path2 <- paste0("data/intermediate/boot_support.rds")
tr <- readRDS(file = input_path1)
bt <- readRDS(file = input_path2)

source("code/functions/bs_analysis.r")
source("code/parameters/boot.r")

#plot bs analisis
for (m in 1:3){
pdf(file = paste0("data/intermediate/bs_analysis_mode", m, ".pdf"),
 height = 8, width = 10)
par(mfrow = c(2, 2))
for (p in 1:length(tr)) bs_analysis(tr[p], bt[[p]], mode = m)
dev.off()
}
