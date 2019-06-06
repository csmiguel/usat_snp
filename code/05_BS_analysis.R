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

bs_treemetrics <-
  seq_along(tr) %>%
    lapply(function(x) tree_metrics(tr[[x]], bt[[x]]))
names(bs_treemetrics) <- names(tr)

saveRDS(bs_treemetrics, "data/intermediate/bs_treemetrics.rds")
  #statistics
  dnnmin.glm <- summary(glm(bs~dnnmin, family = binomial()))
    dnnmin.glm.e <- dnnmin.glm$coefficients[2, 1] %>% round(3)
    dnnmin.glm.s <- dnnmin.glm$coefficients[2, 4] %>%
      formatC(format = "e", digits = 2)
  dno.glm <- summary(glm(bs~dno, family = binomial()))
    dno.glm.e <- dno.glm$coefficients[2, 1] %>% round(3)
    dno.glm.s <- dno.glm$coefficients[2, 4] %>%
      formatC(format = "e", digits = 2)
