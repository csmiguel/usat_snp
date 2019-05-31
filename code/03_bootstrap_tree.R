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
#    Description:
#    Outpath:
#  DEPENDENCIES:
###.............................................................................
library(adegenet)
library(dplyr)
library(ape)

input_path <- "data/intermediate/gen_consolidated_filtered.rds"
gen <- readRDS(file = input_path)

source("code/functions/bootstrap_trees.r")
source("code/functions/create_trees.r")
source("code/parameters/boot.r") #number of bootstrap

#estimate tree from original data and from bootstrapped genotypes
tr <- seq_along(gen) %>%
  sapply(function(j){
    if (adegenet::nLoc(gen[[j]]) < 100) type <- "usat"
    if (adegenet::nLoc(gen[[j]]) > 100) type <- "snp"
    sp <- names(gen)[j]
    temp <- list()
    #create nj tree on original genotypes
    original_trees <- create_trees(gen[[j]], type = type)
    #create trees from bootstrapped genotypes
    boots_trees <- bootstrap_nj(gen[[j]], nboot, type = type)
    temp[[1]] <- c(original_trees, boots_trees)
    names(temp) <- sp
    temp
})

#root trees
tr_rooted <- midpoint_tr(tr)

#get bootstrap
bt <- get_bt(tr_rooted)

#save trees (first tree is the real tree)
saveRDS(tr, "data/intermediate/bt_trees.rds")
saveRDS(bt, "data/intermediate/boot_support.rds") #bootstrap support
