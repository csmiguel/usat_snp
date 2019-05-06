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

input_path1 <- paste0("data/intermediate/filt_genotypes.rds")
gen <- readRDS(file = input_path1)
input_path2 <- paste0("data/intermediate/bt_trees.rds")
tr <- readRDS(file = input_path2)
input_path3 <- paste0("data/intermediate/boot_support.rds")
bt <- readRDS(file = input_path3) %>% sapply(function(x) x / 100)

source("code/parameters/boot.r")
source("code/functions/bootstrap_trees.r")
source("code/functions/create_trees.r")

#subset data
pelo <- gen$dart_pelo
original_trees <- lapply(tr, function(x) x[[1]])
rm(tr)

#create trees from subsampled dataset
s.tr <- 1:length(s) %>%
  sapply(function(j){
  sp <- as.character(s[j])
  temp <- list()
  for (i in 1:n_sampling_replicas){
    int <- sample(1:nLoc(pelo), s[j], replace = FALSE)
    s_pelo <- pelo[, int]
    if (i == 1)
    boots_trees <- bootstrap_nj(s_pelo, nboot, type = "snp")
    if (i > 1)
    boots_trees <- c(boots_trees, bootstrap_nj(s_pelo, nboot, type = "snp"))
  }
  temp[[1]] <- boots_trees
  names(temp) <- sp
  temp
})

#calculate bootstrap support
s.bt <- names(s.tr) %>% seq_along() %>%
  sapply(function(x){
    sp <- names(s.tr[x])
    temp <- list()
    phy <- original_trees$dart_pelo
    ints <- phy$edge[, 2] > ape::Ntip(phy)
    ans <- ape::countBipartitions(phy, s.tr[[x]])
    temp[[1]] <- c(nboot * n_sampling_replicas, ans[order(phy$edge[ints, 2])])
    names(temp) <- sp
    temp
})

#create data frame with combined BS for subsampled trees and original ones
s.btx <- lapply(s.bt, function(x) x[-1])
s_bt_m <- reshape::melt(s.btx)
s_bt_m[, 1] <- s_bt_m[, 1] / (nboot * n_sampling_replicas)
assertthat::assert_that(max(s_bt_m[, 1]) == 1)
s_bt_m[, 2] <- as.numeric(s_bt_m[, 2]) %>% as.factor()
s_bt_m <- rbind(s_bt_m, reshape::melt(bt))
names(s_bt_m) <- c("BS", "dataset")
s_bt_m$BS <- round(s_bt_m$BS, 2)

#quantiles for BS
g <- by(data = s_bt_m$BS, INDICES = s_bt_m$dataset, quantile)
gg <- 1:length(g) %>% sapply(function(x) g[[x]] %>% round(2))
colnames(gg) <- names(g); rm(g)

#labs for plot
nl <- sapply(gen, adegenet::nLoc)
lab <- c(levels(s_bt_m$dataset)[1:length(s)],
  paste0(names(gen), rep("\n (", 4), nl, rep(")", 4)))

#boxplot
pdf(file = "data/intermediate/bs_boxplot.pdf", height = 5, width = 8)
  boxplot(s_bt_m$BS~s_bt_m$dataset, xlab = "Number of loci",
      ylab = "Bootstrap support",
      main = "Change of overall boostrap support accross number of loci",
      col = c(rep("grey", length(s)), rep(grey(0.3), 2), rep(grey(0.8), 2)),
      pch = 20, outcol = grey(0.7), xaxt = "n", yaxt = "n")
      axis(1, at = 1:(length(s) + 4), labels = lab, las = 2)
      axis(2, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), las = 1)
      abline(v = length(s) + 0.5, lty = 2, lwd = 2)
dev.off()


#save trees (first tree is the real tree)
saveRDS(s.tr, "data/intermediate/s_trees.rds")
saveRDS(gg, "data/intermediate/s_quantiles.rds")
