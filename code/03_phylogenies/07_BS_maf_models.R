###.............................................................................
# (c) Miguel Camacho-Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: fit models for relation between BS no loci
#DESCRIPTION:
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

library(dplyr)

raw_data <- readRDS("data/intermediate/maf_bs_treemetrics.rds")

# Preliminary processing
## prepare data

species <- sapply(strsplit(names(raw_data), "_", fixed = TRUE), `[[`, 2)

data <-
  plyr::ldply(seq_along(raw_data), function(i){
    plyr::ldply(seq_along(raw_data[[i]]), function(j){
      x <- raw_data[[i]][[j]]
      x$maf <- names(raw_data[[i]])[j]
      x$species <- species[i]
      x
    })
  }) %>%
    dplyr::select(-dnp, -dnr) %>%
  dplyr::mutate(maf = as.numeric(maf)) %>%
  as_tibble()

# fit models
maf_models <-
  plyr::dlply(data, c("species"), function(x) {
  suppressWarnings({
    full <- lm(bs_tm ~ maf, data = x)
    null <- lm(bs_tm ~ 1, data = x)
  })
  list(full = full, null = null)
})

# test for difference between number of loci
maf_test <-
  lapply(maf_models, function(x) {
  suppressWarnings(broom::tidy(anova(x$full, x$null, test = "Chisq")))
})

#sink results models
sink("data/final/BS_models_MAF.txt")

cat("Results from models:\n1.1. MAF models")
lapply(maf_models, function(x) lapply(x, summary))
cat("\n\n1.2. MAF tests")
maf_test

sink()
