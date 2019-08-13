###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: subsampling snps bs support
#DESCRIPTION:
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

library(magrittr)
library(dplyr)
library(ape)
library(dartR)

#genotypes
input_path1 <- "data/intermediate/shared_all_consolidated_filtered.rds"
gen <- readRDS(file = input_path1) %>%
  {.[grepl("dart", names(.))]} #only select dart datasets
#trees
input_path2 <- "data/intermediate/shared_bt_trees.rds"
tr <- readRDS(file = input_path2)

#tree metrics from datasets
input_path3 <- "data/intermediate/shared_bs_treemetrics.rds"
bs <- readRDS(input_path3)

source("code/parameters/boot.r")
source("code/functions/bootstrap_trees.r")
source("code/functions/bs_analysis.r")

#empty list
out_trees <- list()

#calculate bootstrap support with tree metrics (dnr and dnp) for trees made from
#datasets which are a subsample of the loci from larger matrices.
bs_treemetrics.s <-
seq_along(gen) %>%
  lapply(function(z){
#   subset data
    gen1 <- gen[[z]] %>% dartR::gi2gl()
    original_tree <- lapply(tr, function(x) x[[1]])[[names(gen)[z]]]

#   create trees from subsampled dataset
#   no of bootstrap per loci sampling
    nboots <- nboot / n_sampling_replicas
#   limit the number of subsampled loci to the a number below total loci
    ss <- s[s < adegenet::nLoc(gen1)]
#   create bootstraped trees
      s.tr <-
      seq_along(ss) %>%
      sapply(function(j){
        sp <- as.character(ss[j])
        temp <- list()
        #for each loci subsampling
        for (i in 1:n_sampling_replicas){
          #subsample s[j] loci without replacement
          int <- sample(1:adegenet::nLoc(gen1), ss[j], replace = FALSE)
          s_gen1 <- gen1[, int]
          #create nboots replicate trees
          if (i == 1)
           boots_trees <- bootstrap_nj(s_gen1, nboots, type = "snp")
          if (i > 1)
           boots_trees <- c(boots_trees,
             bootstrap_nj(s_gen1, nboots, type = "snp"))
        }
        temp[[1]] <- boots_trees
        names(temp) <- sp
        temp
      })
  out_trees[[z]] <<- s.tr
  #root trees
  tr_rooted <- midpoint_tr(s.tr)

  #get bootstrap#the function internally assigns 0's to all NA values in
  #bootstrap calculations.
  bt <- get_bt(tr_rooted, external_ref_tree = original_tree)
  h <-
  seq_along(tr_rooted) %>%
    lapply(function(x){
      tree_metrics(tr_rooted[[x]], bt[[x]], external_ref_tree = original_tree)
          })
      names(h) <- names(s.tr)
      h
      })
names(bs_treemetrics.s) <- names(gen)
names(out_trees) <- names(gen)

saveRDS(out_trees, "data/intermediate/shared_s_trees.rds")
saveRDS(bs_treemetrics.s, "data/intermediate/shared_s_bs_treemetrics.rds")
