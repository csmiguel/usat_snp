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

library(dplyr)
library(ape)
library(dartR)

#genotypes
input_path <- "data/intermediate/shared_all_consolidated_filtered.rds"
gen <- readRDS(file = input_path) %>%
  {.[grepl("dart", names(.))]} #only select dart datasets

#trees
input_path2 <- "data/intermediate/shared_bt_trees.rds"
tr <- readRDS(file = input_path2)

source("code/parameters/boot.r")
source("code/functions/bootstrap_trees.r")
source("code/functions/bs_analysis.r")

#empty list
out_trees <- list()

#calculate bootstrap support across different mafs

bs_treemetrics.maf <-
  seq_along(gen) %>%
  lapply(function(z){

    original_tree <- lapply(tr, function(x) x[[1]])[[names(gen)[z]]]

#   convert from genind to genlight
    gen1 <- gen[[z]]
#   remove objects within
    other(gen1) <- NULL
#   regenerate $loc.metrics to be read correctly by dartR functions.
    gen1 %<>% dartR::gi2gl() %>%
      dartR::gl.recalc.metrics()

#   filter pelobates genotypes accross MAF thresholds
    gen_maf <-
      seq_along(maf) %>%
      lapply(function(x){
        dartR::gl.filter.maf(gen1, threshold = maf[x], v = 5)
      }
      ) %>% {
        c(gen1, .)#concatenate original genotypes (maf = 0.02) to filtered
      }
#    vector with all maf thresholds
    all_maf <- as.character(c(0.02, maf))
    names(gen_maf) <- all_maf

#   create trees from maf datasets
    nboots <- nboot / n_sampling_replicas

    maf.tr <-
      seq_along(all_maf) %>%
      sapply(function(j){
        sp <- all_maf[j]
        temp <- list()
        pmaf <- gen_maf[[j]]
        for (i in 1:n_sampling_replicas){
          int <- sample(1:adegenet::nLoc(pmaf), mafss, replace = FALSE)
          s_pmaf <- pmaf[, int]
          if (i == 1)
            boots_trees <- bootstrap_nj(s_pmaf, nboots, type = "snp")
          if (i > 1)
            boots_trees <- c(boots_trees,
              bootstrap_nj(s_pmaf, nboots, type = "snp"))
        }
        temp[[1]] <- boots_trees
        names(temp) <- sp
        temp
      })
    out_trees[[z]] <<- maf.tr

#   root trees
    tr_rooted <- midpoint_tr(maf.tr)

#   get bootstrap#the function internally assigns 0's to all NA values in
#   bootstrap calculations.
    bt <- get_bt(tr_rooted, external_ref_tree = original_tree)

#   tree metrics
    h <-
      seq_along(tr_rooted) %>%
      lapply(function(x) tree_metrics(tr_rooted[[x]], bt[[x]],
                                      external_ref_tree = original_tree))
    names(h) <- names(maf.tr)
    h
      })
names(bs_treemetrics.maf) <- names(gen)
names(out_trees) <- names(gen)

#save objecs
saveRDS(out_trees, "data/intermediate/shared_maf_trees.rds")
saveRDS(bs_treemetrics.maf, "data/intermediate/shared_maf_bs_treemetrics.rds")
