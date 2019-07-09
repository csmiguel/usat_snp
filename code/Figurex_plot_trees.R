###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: Plot results from bootstrap analysis
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(ape)
library(dplyr)
library(adegenet)

#bootstrap support
bt <- readRDS("data/intermediate/boot_support.rds")
#trees
tr <- readRDS("data/intermediate/bt_trees.rds")
#genotypes
gen <- readRDS("data/intermediate/gen_consolidated_filtered.rds")
#best files as determined by KFinder
best_files <- readRDS("data/intermediate/best_files_kfinder.rds")
#clumpak output
qclumpak <- readRDS("data/intermediate/clumpak_major.rds")
#pop ids
meta <- readRDS("data/intermediate/table1.rds")

source("code/parameters/boot.r")
source("code/functions/plot_nj_trees.r")
#rename datasets
names(tr) <-
  plyr::mapvalues(
    x = names(tr),
    from = names(tr),
    to = c("SNPs H. molleri", "SNPs P. cultripes",
    "microtellites H. molleri", "microtellites P. cultripes"))

#plot trees
plot_nj(bt = bt, thresh = thresh, coln = c("red", "black"),
  tips = "none", typeP = "unrooted", K = "K2", output = ".pdf")
#plots using ancestries from major cluster in K = 2
plot_nj(bt = bt, thresh = thresh, coln = c("red", "black"),
  tips = "ancestry", typeP = "phylogram", K = "K2", output = ".pdf")
#plots using ancestries from major cluster in K = 4
plot_nj(bt = bt, thresh = thresh, coln = c("red", "black"),
  tips = "ancestry", typeP = "phylogram", K = "K4", output = ".pdf")
