###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: compute standarized multilocus heterozygosity
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(dartR)
library(dplyr)

#genotypes
gen <- readRDS("data/intermediate/gen_consolidated_filtered.rds")

source("code/functions/genind2inbreedR.r") #genind2 inbreedR

##
het <-
  lapply(names(gen), function(z){
  #compute and map heterozygosity
  l <- dartR::gi2gl(gen[[z]]) %>% #select genetic dataset
       as.matrix()
  l[l == 2] <- 0 #reformat for inbreedR: homozygotes == 0.
  if (grepl("usat", z)) #for microsatellites reformat matrix to compute ht.
    l <- usatgenind2inbreedR(gen[[z]])
  assertthat::assert_that(check_data2(l))
  #standarized heterozygosity
  het2 <- inbreedR::sMLH(l)
})
names(het) <- names(gen)

saveRDS(het, "data/intermediate/sMLH.rds")
