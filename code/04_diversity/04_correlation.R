###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: compute correlation between sMLH usats vs SNPs
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
  dplyr::as_tibble()


#pearson correlation
p1.test <-
  plyr::dlply(data, ~species, function(x){
    cor.test(x$usats, x$snps, method = "pearson")
    })

#write results from correlation
sink("data/final/correlation_sMLH.txt")
p1.test
sink()
