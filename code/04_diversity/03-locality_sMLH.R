###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: median sMLH per locality
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(dartR)
library(dplyr)

#read sMLH
het <- readRDS("data/intermediate/sMLH.rds")
#genotypes
gen <- readRDS("data/intermediate/gen_consolidated_filtered.rds")

#assertion
assertthat::assert_that(all(names(het) %in% names(gen)))

#compute median sMLH per locality

median_het <-
  names(het) %>%
    lapply(function(x){
      assertthat::assert_that(all(names(het) %in% names(gen)))
      #sort heterozygosity by names in
      h <- het[[x]][gen[[x]]$other$metadata$sample_id]
      split(h, gen[[x]]$other$metadata$locality) %>%
        sapply(median)
      }) %>% setNames(names(het))

saveRDS(median_het, "data/intermediate/median_het.rds")
