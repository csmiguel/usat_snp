###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: Table with samples
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(dplyr)
library(flextable)

#load consolidated genotypes
gen <- readRDS("data/intermediate/shared_gen_consolidated_filtered.rds")
#load all gen
gen_all <- readRDS("data/intermediate/gen_consolidated_filtered.rds")
#gen with shared for Pelobates and originals for Hyla
seq_along(gen) %>%
  lapply(function(x){
    gen_all[[match(names(gen)[x], names(gen_all))]] <<- gen[[x]]
    })
gen <- gen_all
rm(gen_all)

source("code/functions/centroid_sf.r")
source("code/functions/ft_helper.r")

#one table with all data
h <-
seq_along(gen) %>%
  lapply(function(x){
    gen[[x]]@other$metadata %>%
      sfc_as_cols %>% #sf object to dataframe with lat lon columns
      dplyr::mutate(dataset = names(gen)[x]) %>% #add a column with dataset
      #round coordinates
      dplyr::mutate(longitude = round(longitude, 3)) %>%
      dplyr::mutate(latitude = round(latitude, 3))
  }) %>%
  plyr::ldply(rbind) %>% #bind all tables within the list
  tibble::as_tibble()

#inspect duplicated samples
dd <- h$sample_id[which(duplicated(h$sample_id))] #dublicated ids
#arrange by sample id
h %>% dplyr::filter(sample_id %in% dd) %>% dplyr::arrange(sample_id)

#there are some incongruences in the coordinates (difference in a hundreth of a degree) for some of the duplicated samples. I have checked at the coordinates on the original metadata and the publications and everything seems to be fine, so I do not where the error comes from.

#table 1: expanded including duplicates
p <-
  h %>%
  dplyr::mutate(species = ifelse(grepl("hyla", dataset), "H.m.",
                                ifelse(grepl("pelo", dataset), "P.c.", NA))) %>%
  dplyr::mutate(marker = ifelse(grepl("dart", dataset), "SNPs",
                     ifelse(grepl("usat", dataset), "microsatellites", NA))) %>%
  select(-dataset) %>%
  arrange(locality, sample_id)

#table 2: compact

#calculate centroid for each population
localities <-
  p %>%
  #convert to spatial points object
  sp::SpatialPointsDataFrame(coords = cbind(p$latitude, p$longitude)) %>%
  as("sf") %>% #convert to sf
  sf::st_set_crs(4326) %>%
  centroid_sf() %>%
  sfc_as_cols %>%
  dplyr::mutate(longitude = round(longitude, 3)) %>%
  dplyr::mutate(latitude = round(latitude, 3)) %>% tibble::as_tibble()

#factors to iterate
factorss <- expand.grid(species = unique(p$species), marker = unique(p$marker))
#create table 2 with uniques rows per locality and sample ids
pp <-
1:nrow(factorss) %>%
  lapply(function(y){
    sp <- factorss[y, 1]
    mk <- factorss[y, 2]
    p %>%
      dplyr::filter(marker == mk & species == sp) %>%
      dplyr::arrange(locality, species, marker) %>%
      dplyr::select(-longitude, -latitude) %>%
      aggregate(by = list(.$locality), function(x) x) %>%
      dplyr::select(1, sample_id) %>%
      dplyr::mutate(sample_id =
                  sapply(sample_id, function(x) paste(x, collapse = ", "))) %>%
      mutate(sample_id = paste0(sp, ", ", mk, ": ", sample_id))
    })

#final formatting
ppp <-
  localities %>%
  dplyr::left_join(pp[[1]], by = c("locality" = "Group.1")) %>%
  dplyr::left_join(pp[[2]], by = c("locality" = "Group.1")) %>%
  dplyr::left_join(pp[[3]], by = c("locality" = "Group.1")) %>%
  dplyr::left_join(pp[[4]], by = c("locality" = "Group.1")) %>%
  #paste samples from all columns
  tidyr::unite(new, sample_id.x, sample_id.y,
               sample_id.x.x, sample_id.y.y, sep = "; ") %>%
  #replacements
  dplyr::mutate(samples = gsub(pattern = "NA; |; NA",
    replacement = "", x = new)) %>%
  dplyr::select(-new) %>%
  #create unique IDs per population
  tibble::rowid_to_column("ID")
#export as word
ft <-
ppp %>% flextable::flextable() %>%
  flextable::width(j = 5, 3.7) %>%
  flextable::width(j = 2, 2) %>%
  flextable::fontsize(size = 8, part = "body") %>%
  flextable::fontsize(size = 10, part = "header")
ft

ft2word(ft = ft, name = "data/final/shared_table1.docx", ftcaption = c("Table 1. Samples and their localities. In the sample column, we indicate the species(*P.c.* and *H.m*, to refer to *Pelobates cultripes* and *Hyla molleri*, respectively) and marker type (SNPs or microsatellite)."))

saveRDS(ppp, "data/intermediate/shared_table1.rds")
saveRDS(gen, "data/intermediate/shared_all_consolidated_filtered.rds")
