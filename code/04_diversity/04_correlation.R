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
  dplyr::select(species, locality, sMLH_SNPs_median, sMLH_usats_median) %>%
  dplyr::rename(snps = sMLH_SNPs_median,
                usats = sMLH_usats_median) %>%
  dplyr::mutate(species = as.factor(species),
                locality = as.factor(locality)) %>%
  dplyr::as_tibble()


#pearson correlation
p1.test <-
  plyr::dlply(data, ~species, function(x){
    cor.test(x$usats, x$snps, method = "pearson")
    })

#linear model
m1 <-
  plyr::dlply(data, ~species, function(x){
    lm(x$snps ~ x$usats)
  })

#residuals
#get residuals from models
res <-
  seq_along(m1) %>%
  lapply(function(x){
    m1[[x]]$residuals #residuals from 1 sp
    }) %>% setNames(levels(data$species))

#for each species attach residuals to data
h <- plyr::dlply(data, ~species)
resAll <-
  seq_along(h) %>%
  lapply(function(x){
    hh <- h[[x]]
    hh$res <- round(res[[x]], 2)
    hh
  }) %>% setNames(levels(data$species))

#for P. cultripes compute median per population
resPelo <-
  split(resAll$`P. cultripes`$res, resAll$`P. cultripes`$locality, drop = T) %>%
  sapply(median) %>%
  {data.frame(res = .)} %>%
  tibble::rownames_to_column("locality") %>%
  dplyr::mutate(species = levels(data$species)[2])

#merge to Hyla (already median per locality)
residuals_formatted <-
  resAll$`H. molleri` %>%
  dplyr::select(locality, res, species) %>%
  rbind(resPelo) %>%
  dplyr::select(3, 1, 2)

#write results from correlation
sink("data/final/correlation_sMLH.txt")
p1.test
sink()

#save residuals from linear model
saveRDS(residuals_formatted, "data/intermediate/residuals-sMLH.rds")
