###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: make marker datasets for each species more comparable:
#      1. Hyla: keep only samples from close localities ,~< 30 km
#      2. Pelobates: retain samples shared between datasets.
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

#load libraries
library(dartR)
library(dplyr)
library(magrittr)

#load data and functions
input_path <- "data/intermediate/gen_consolidated.rds"
gen <- readRDS(file = input_path)
source("code/functions/map_distribution.r")

###########
# 1. Hyla #
###########

#Plot samples in map so as to decide which populations to keep
plot_samples(meta = gen,
  name_plot = "data/intermediate/pre_sample_fil_distribution.pdf", res = 5000)

#populations to remove based on map:
rm_pop <-
  list(
    usat_hyla =
      c("Segonnes", "Larués", "Renales", "Torrecaballeros",
        "Valdemanco", "Colmenar Viejo", "Peñalara", "Río Montoro",
        "Albires", "Serra de Montejunto", "Mindelo", "El Puerto", "Valdeazores",
        "Valle de Santiago", "Boticas", "Fuenterrebollo", "Ciudad Rodrigo"),
    dart_hyla = NULL)

#filter localities
assertthat::assert_that(all(names(rm_pop) %in% names(gen)))
for (i in names(rm_pop)){
  gen[[i]] <- gen[[i]][!pop(gen[[i]]) %in% rm_pop[[i]]]
  assertthat::assert_that(
    nrow(gen[[i]]$other$metadata) == adegenet::nInd(gen[[i]]))
}
rm(i)

#assert all of the pop names == metadata$locality names
gen %>% sapply(function(x){
  assertthat::assert_that(
    all(
      x$other$metadata$locality == pop(x)
      ))})

#edit locality names in metadata to match localities between maker datasets
# Hyla:
# 1. change names
gen$dart_hyla$other$metadata$locality %<>%
  gsub(pattern = "Xinzo de Limia", replacement = "Ginzo de Limia") %>%
  gsub(pattern = "Laguna Ojos de Villaverde",
    replacement = "Ojos de Villaverde") %>%
  gsub(pattern = "Arzila", replacement = "Lavariz") %>%
  gsub(pattern = "Serra da Estrella", replacement = "Serra da Estrela") %>%
  gsub(pattern = "Valgañon", replacement = "Valgañón") %>%
  gsub(pattern = "Casa Velha", replacement = "Longueira") %>%
  gsub(pattern = "Sao Jose de Lamarosa|Coruche", replacement = "Santarém")
gen$usat_hyla$other$metadata$locality <-
  gsub("Castelo de Vide", "Beira", gen$usat_hyla$other$metadata$locality)
# 2. merge localities
gen$dart_hyla$other$metadata$locality %<>%
  gsub(pattern = "Saceruela|Río Esteras",
    replacement = "Sierra de Saceruela")
gen$usat_hyla$other$metadata$locality %<>%
  gsub(pattern = "Fontanosas", replacement = "Sierra de Saceruela")
# 2.2
gen$usat_hyla$other$metadata$locality %<>%
  gsub(pattern = "Los Baños de Robledillo|Navas de Estena",
    replacement = "Navas de Estena - Baños Robledillo")
gen$dart_hyla$other$metadata$locality %<>%
  gsub(pattern = "Navas de Estena",
    replacement = "Navas de Estena - Baños Robledillo")
# 2.3
gen$dart_hyla$other$metadata$locality %<>%
  gsub(pattern = "Codesal|Ferreras de Arribas",
    replacement = "Codesal - Ferreras de Arriba")
gen$usat_hyla$other$metadata$locality %<>%
  gsub(pattern = "Codesal|Ferreras de Arriba",
    replacement = "Codesal - Ferreras de Arriba")

#Pelobates:
# 1. change names
gen$dart_pelo$other$metadata$locality %<>%
  gsub(pattern = "Plage Biscarrose", replacement = "Plage Biscarrosse") %>%
  gsub(pattern = "La Porroquia", replacement = "La Parroquia")
# 2. merge localities
gen$usat_pelo$other$metadata$locality %<>%
  gsub(pattern = "Benalup-Casas Viejas", replacement = "Los Alcornocales")
gen$dart_pelo$other$metadata$locality %<>%
  gsub(pattern = "Tarifa", replacement = "Los Alcornocales")


################
# 2. Pelobates #
################
#Keep only shared individuals for Pelobates

# vector with individuals to keep
samples_shared <-
  adegenet::indNames(gen$dart_pelo)[(
    adegenet::indNames(gen$dart_pelo) %in% adegenet::indNames(gen$usat_pelo))]

# filter shared samples
for (i in grep("pelo", names(gen), value = T)){
  gen[[i]] <-
    gen[[i]][match(samples_shared, adegenet::indNames(gen[[i]]))]
}
rm(i)

# check order and number of individuals between markers are identical
assertthat::assert_that(all(
  adegenet::indNames(gen[["usat_pelo"]]) ==
  adegenet::indNames(gen[["dart_pelo"]])))


###################################
# Finalize metadata consolidation #
###################################

#consolite pop slot
gen %<>% lapply(function(x){
  #assert metadata$sample_id are the same to indNames
  assertthat::assert_that(
    all(x$other$metadata$sample_id == adegenet::indNames(x)))
  pop(x) <- x$other$metadata$locality
  assertthat::assert_that(
    identical(x$other$metadata$locality, as.character(pop(x))))
  x
  })

#plot localities of filtered datasets
plot_samples(meta = gen,
  name_plot = "data/intermediate/post_sample_fil_distribution.pdf", res = 5000)

# assert localities in dartseq and usat are the same
assertthat::assert_that(
  #Hyla
  identical(
  levels(as.factor(gen$usat_hyla$other$metadata$locality)),
  levels(as.factor(gen$dart_hyla$other$metadata$locality))) &
  #Pelobates
  identical(
    levels(as.factor(gen$usat_pelo$other$metadata$locality)),
    levels(as.factor(gen$dart_pelo$other$metadata$locality)))
  )

#save object
saveRDS(gen, "data/intermediate/gen_consolidated_filtered.rds")
