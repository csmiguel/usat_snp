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
library(dplyr)
library(ape)

input_path1 <- "data/intermediate/bt_trees.rds"
input_path2 <- "data/intermediate/boot_support.rds"
tr <- readRDS(file = input_path1)
bt <- readRDS(file = input_path2)

source("code/functions/bs_analysis.r")
source("code/parameters/boot.r")

#get list of dataframes with columns being:
# dnr, distance from node to root
# dnp, distance from node to parent node
# bs_tm, bootstrap support
#and rows being nodes.

bs_treemetrics <-
  seq_along(tr) %>%
    lapply(function(x) tree_metrics(tr[[x]], bt[[x]]))
names(bs_treemetrics) <- names(tr)

saveRDS(bs_treemetrics, "data/intermediate/bs_treemetrics.rds")
