###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: compute residuals from regression sMLH snps vs usats with slope = 1.
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(dplyr)

#load sMLH
het <- readRDS("data/intermediate/sMLH_reformatted.rds")

#reformat data. For Hyla keep only localities and for Pelobates keep all ind.
data <-
  het %>%
  #remove rows with duplicated localities for Hyla
  {.[-which(.$species == "H. molleri" & duplicated(.$locality)), ]} %>%
  dplyr::mutate(snps = ifelse(species == "H. molleri", sMLH_SNPs_median,
                        sMLH_SNPs),
                usats = ifelse(species == "H. molleri", sMLH_usats_median,
                                      sMLH_usats)) %>%
  dplyr::select(species, locality, snps, usats) %>%
  dplyr::mutate(species = as.factor(species),
                locality = as.factor(locality)) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(res = snps - usats)

saveRDS(data, "data/intermediate/sMLH_residuals.rds")
