###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: analysis of bootstrap support for (1) all data, and for (2) comparison
# between usats and snps. For goal (2) I will predict the number of SNPs that
# are needed for each species that is needed to meet the average BS for
# usats data for that given species.
#DESCRIPTION:
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
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

# (1) all data
#get list of dataframes with columns being:
# dnr, distance from node to root
# dnp, distance from node to parent node
# bs_tm, bootstrap support
#and rows being nodes.

bs_treemetrics <-
  seq_along(tr) %>%
    lapply(function(x) tree_metrics(tr[[x]], bt[[x]]))
names(bs_treemetrics) <- names(tr)

# (2) comparison snps-usats
# 2.1 get number of SNPs which have around the same amounts of information than
# usat data for that given species.
subsampling_models <- readRDS("data/intermediate/subsampling_models.rds")

predloci <-
  seq_along(subsampling_models) %>%
  sapply(function(sp){
    btusat <- bt[[paste0("usat_", names(subsampling_models)[sp])]][-1] %>%
      mean / 1000 #average bs usat
    coe <- subsampling_models[[sp]]$full$coefficients #coefficients of the model
    #compute nloci needed to get a BS support in the SNP data similar to usats
    exp( (btusat - coe["(Intercept)"]) / coe["log(nloci)"]) %>% round()
    }) %>% setNames(paste0("dart_", names(subsampling_models)))

predloci80 <-
  seq_along(subsampling_models) %>%
  sapply(function(sp){
    bt80 <- .8  #average bs to predict
    coe <- subsampling_models[[sp]]$full$coefficients #coefficients of the model
    #compute nloci needed to get a BS support in the SNP data similar to usats
    exp( (bt80 - coe["(Intercept)"]) / coe["log(nloci)"]) %>% round()
    }) %>% setNames(paste0("dart_", names(subsampling_models)))

# 2.2 get metrics

#load data genotypes
input_path1 <- "data/intermediate/gen_consolidated_filtered.rds"
gen <- readRDS(file = input_path1) %>%
  {.[grepl("dart", names(.))]} #only select dart datasets

# load functions
source("code/parameters/boot.r")
source("code/functions/bootstrap_trees.r")
source("code/functions/bs_analysis.r")
rm(s)

#empty list
out_trees <- list()

#calculate bootstrap support with tree metrics (dnr and dnp) for trees made from
#datasets which are a subsample of the loci from larger matrices.
bs_treemetrics.s <-
seq_along(predloci) %>%
  lapply(function(z){
#   subset data
    gen1 <- gen[[names(predloci)[z]]]
    original_tree <- lapply(tr, function(x) x[[1]])[[names(gen)[z]]]

#   create trees from subsampled dataset
#   no of bootstrap per loci sampling
    nboots <- nboot / n_sampling_replicas
#   limit the number of subsampled loci to the a number below total loci
    ss <- predloci[names(predloci)[z]]
#   create bootstraped trees
      s.tr <-
      seq_along(ss) %>%
      sapply(function(j){
        sp <- as.character(ss[j])
        temp <- list()
        #for each loci subsampling
        for (i in 1:n_sampling_replicas){
          #subsample s[j] loci without replacement
          s_gen1 <- gen1[loc = sample(locNames(gen1), ss[j], replace = F)]
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

# merge data from usats to subsamples dart
bs_treemetrics_comp <-
  lapply(bs_treemetrics.s, function(x) x[[1]]) %>%
  c(bs_treemetrics[grep("usat", names(bs_treemetrics))])

# save objects
# all data
saveRDS(bs_treemetrics, "data/intermediate/bs_treemetrics.rds")
# data for comparison
saveRDS(bs_treemetrics_comp, "data/intermediate/bs_treemetrics_comp.rds")
saveRDS(predloci, "data/intermediate/predloci.rds")
saveRDS(predloci80, "data/intermediate/predloci.rds")
