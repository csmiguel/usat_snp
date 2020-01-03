###.............................................................................
# (c) Jeffrey Hanson developped most of this script
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: fit models for relation between BS and markers (for each sp and metric)
#DESCRIPTION:
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

library(dplyr)
library(tibble)
library(tidyr)

# 1. Compare BS averages between makers for each species using all data

#read dnp, dnr and bootstrap support from subset of comparable SNPs
raw_data_all <- readRDS("data/intermediate/bs_treemetrics.rds")

# Preliminary processing
## prepare data
marker <- sapply(strsplit(names(raw_data_all), "_", fixed = TRUE), `[[`, 1)
species <- sapply(strsplit(names(raw_data_all), "_", fixed = TRUE), `[[`, 2)
data_all <-
  plyr::ldply(seq_along(raw_data_all), function(i) {
    x <- raw_data_all[[i]]
    x$marker <- marker[i]
    x$species <- species[i]
    x
  }) %>%
  as_tibble()

bs_models <- plyr::dlply(data_all, c("species"), function(x) {
  suppressWarnings({
    glm(bs_tm ~ marker, data = x, family = "binomial") %>% summary()
  })
})

sink("data/final/BS_models_means_marker_all.txt")
bs_models
sink()

rm(list = ls())
# 2. Compare change of Dnp and Dnr between markers for each species

#read dnp, dnr and bootstrap support from subset of comparable SNPs
raw_data <- readRDS("data/intermediate/bs_treemetrics_comp.rds")

# Preliminary processing
## prepare data
marker <- sapply(strsplit(names(raw_data), "_", fixed = TRUE), `[[`, 1)
species <- sapply(strsplit(names(raw_data), "_", fixed = TRUE), `[[`, 2)
data <- plyr::ldply(seq_along(raw_data), function(i) {
  x <- raw_data[[i]]
  x$marker <- marker[i]
  x$species <- species[i]
  x
}) %>%
  as_tibble()

## fit models
models <- plyr::dlply(data, c("species"), function(x) {
  suppressWarnings({
    mbinomial <- glm(bs_tm ~ dnp * marker + dnr * marker, data = x, family = "binomial")
    mgaussian <- glm(bs_tm ~ dnp * marker + dnr * marker, data = x, family = "gaussian")
  })
  list(binomial = mbinomial, gaussian = mgaussian)
})
#compare models
glm.diag.plots(m1)
glm.diag.plots(m2)

## format raw data for plotting
plot_data <-
  data %>%
  tidyr::gather(dist, value, -bs_tm, -marker, -species)

#export objects
saveRDS(plot_data, "data/intermediate/model_marker_pred_data.rds")

#sink results models
sink("data/final/BS_models_dnpr_marker.txt")
cat("Results from models:\n")
lapply(models, function(x) lapply(x, summary))
sink()
