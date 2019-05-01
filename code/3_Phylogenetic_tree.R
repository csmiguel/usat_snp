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
source("code/functions/create_trees.r")

nboot <- 100 #number of bootstrap

tr <- names(gen) %>% seq_along() %>%
sapply(function(j){
  if (adegenet::nLoc(gen[[j]]) < 100) type <- "usat"
    else if (adegenet::nLoc(gen[[j]]) > 100) type <- "snp"
  sp <- names(gen[j])
  temp <- list()
  original_trees <- create_trees(gen[[j]], type = type)
  boots_trees <- bootstrap_nj(gen[[j]], nboot, type = type)
  temp[[1]] <- c(original_trees, boots_trees)
  names(temp) <- sp
  temp
})


#calculate bootstrap support
bt <- names(tr) %>% seq_along() %>%
sapply(function(x){
  sp <- names(tr[x])
  temp <- list()
  temp[[1]] <- ape::countBipartitions(tr[[x]][[1]], tr[[x]][-1])
  names(temp) <- sp
  temp
})

#save trees (first tree is the real tree)
saveRDS(tr, "data/intermediate/bt_trees.rds")
saveRDS(bt, "data/intermediate/boot_support.rds") #bootstrap support
