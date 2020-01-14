###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: reformat sMLH in one dataframe
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

library(dplyr)
library(dartR)
library(magrittr)

median_het <- readRDS("data/intermediate/median_het.rds") %>%
  lapply(function(x){
    data.frame(median_het = x) %>%
      tibble::rownames_to_column(var = "locality")
  })

het <- readRDS("data/intermediate/sMLH.rds") %>%
  lapply(function(x){
  data.frame(het = x) %>%
    tibble::rownames_to_column(var = "sample_id")
    })

meta <- readRDS("data/intermediate/gen_consolidated_filtered.rds") %>%
  lapply(function(x){
    x@other$metadata %>%
      dplyr::as_tibble() %>%
    dplyr::select(-geometry)
    })

diversity <-
  names(meta) %>%
  lapply(function(x){
   dplyr::left_join(meta[[x]], median_het[[x]],
     by = c("locality" = "locality")) %>%
      dplyr::left_join(het[[x]], by = c("sample_id" = "sample_id")) %>%
      dplyr::rename(sMLH = het) %>%
      dplyr::rename(median_sMLH = median_het)
  }) %>% setNames(names(meta))

species <- c("hyla", "pelo")

diversity_reformatted <-
  species %>%
  lapply(function(sp){ #for each species
  samples <- #dataframe with sample ids for both markers
    diversity[grep(sp, names(diversity))] %>%
      {h <<- do.call(rbind, .)} %>% #bind df
      {.$sample_id} %>%
      unique %>%
      data.frame(sample_id = ., stringsAsFactors = F)

  h %<>%#select non-duplicated samples with their localities
    {.[!duplicated(.$sample_id), ]} %>%
    dplyr::select(sample_id, locality)

  dplyr::left_join(samples, h, by = "sample_id") %>%
  left_join(diversity[[paste0("dart_", sp)]],
            by = c("sample_id"), suffix = c("", "_snps")) %>%
  left_join(diversity[[paste0("usat_", sp)]],
            by = "sample_id", suffix = c("", "_usats")) %>%
    dplyr::mutate(species = plyr::mapvalues(sp, species,
                              c("H. molleri", "P. cultripes"))) %>%
    dplyr::select(species, sample_id, locality, sMLH, median_sMLH, sMLH_usats,
                  median_sMLH_usats) %>%
    dplyr::rename(sMLH_SNPs = sMLH,
                  sMLH_SNPs_median = median_sMLH,
                  sMLH_usats_median = median_sMLH_usats) %>%
    dplyr::arrange(locality) %>%
    dplyr::mutate(sMLH_SNPs = round(sMLH_SNPs, 2),
                  sMLH_usats = round(sMLH_usats, 2),
                  sMLH_SNPs_median = round(sMLH_SNPs_median, 2),
                  sMLH_usats_median = round(sMLH_usats_median, 2))
  }) %>% setNames(species) %>%
  do.call(what = rbind) %>% #combine species
  #fill in gaps from inviduals without median values per population
  plyr::ddply(c("species", "locality"), function(x){
    medsnps <- max(x$sMLH_SNPs_median, na.rm = T)
    x$sMLH_SNPs_median <- medsnps
    medusats <- max(x$sMLH_usats_median, na.rm = T)
    x$sMLH_usats_median <- medusats
    x
  })

write.table(diversity_reformatted,
            file = "data/final/sMLH.txt", quote = F, row.names = F)
saveRDS(diversity_reformatted, "data/intermediate/sMLH_reformatted.rds")
