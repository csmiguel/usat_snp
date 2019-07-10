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

library(magrittr)
library(dplyr)
library(ape)
library(dartR)

#genotypes
input_path1 <- "data/intermediate/gen_consolidated_filtered.rds"
gen <- readRDS(file = input_path1)
#trees
input_path2 <- "data/intermediate/bt_trees.rds"
tr <- readRDS(file = input_path2)

#tree metrics from datasets
input_path3 <- "data/intermediate/bs_treemetrics.rds"
bs <- readRDS(input_path3)

source("code/parameters/boot.r")
source("code/functions/bootstrap_trees.r")
source("code/functions/bs_analysis.r")

#subset data
pelo <- gen$dart_pelo %>% dartR::gi2gl()
original_tree <- lapply(tr, function(x) x[[1]])$dart_pelo
rm(tr)

#create trees from subsampled dataset
nboots <- nboot / n_sampling_replicas

s.tr <-
  seq_along(s) %>%
  sapply(function(j){
    sp <- as.character(s[j])
    temp <- list()
    for (i in 1:n_sampling_replicas){
      int <- sample(1:nLoc(pelo), s[j], replace = FALSE)
      s_pelo <- pelo[, int]
      if (i == 1)
        boots_trees <- bootstrap_nj(s_pelo, nboots, type = "snp")
      if (i > 1)
        boots_trees <- c(boots_trees, bootstrap_nj(s_pelo, nboots, type = "snp"))
    }
    temp[[1]] <- boots_trees
    names(temp) <- sp
    temp
  })

#root trees
tr_rooted <- midpoint_tr(s.tr)

#get bootstrap#the function internally assigns 0's to all NA values in
#bootstrap calculations.
bt <- get_bt(tr_rooted, external_ref_tree = original_tree)

bs_treemetrics.s <-
  seq_along(tr_rooted) %>%
  lapply(function(x) tree_metrics(tr[[x]], bt[[x]],
    external_ref_tree = original_tree))
names(bs_treemetrics.s) <- names(tr_rooted)

saveRDS(s.tr, "data/intermediate/s_trees.rds")
saveRDS(bs_treemetrics.s, "data/intermediate/s_bs_treemetrics.rds")
