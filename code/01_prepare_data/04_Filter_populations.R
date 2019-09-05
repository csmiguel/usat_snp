###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: populations within genotypes to as the usat and SNP datasets to match
#the most.
#DESCRIPTION:
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
#  REQUIRED FILES:
#   Description: genlight objects with raw genotypes for Hyla and Pelobates
#   Inpath: data/intermediate/raw_genotypes.rds
#  OUTPUT:
#    Description: raw genotypes in genlight format
#    Outpath: data/intermediate/filtered_genotypes.rds
###.............................................................................
library(dartR)
library(dplyr)
library(magrittr)

input_path <- "data/intermediate/gen_consolidated.rds"
gen <- readRDS(file = input_path)
source("code/functions/map_distribution.r")

#Plot samples in map so as to decide which populations to keep
plot_samples(meta = gen,
  name_plot = "data/intermediate/pre_sample_fil_distribution.pdf", res = 5000)
#Create list with populations to remove
#based on plotting samples on a map, I removed the following populations
usat_hyla_remove_pop <- c("Segonnes", "Larués", "Renales", "Torrecaballeros",
  "Valdemanco", "Colmenar Viejo", "Peñalara", "Río Montoro",
  "Albires", "Serra de Montejunto", "Mindelo", "El Puerto", "Valdeazores",
  "Valle de Santiago", "Boticas", "Fuenterrebollo", "Ciudad Rodrigo")
dart_hyla_remove_pop <- c() #none filtered
usat_pelo_remove_pop <- c("Valflaunes", "Valdunciel", "Porto Covo",
  "Cabezarrubias del Puerto")
dart_pelo_remove_pop <- c("Pataias", "Gaitán", "Valcárcel")
rm_pop <- paste0(names(gen), "_remove_pop") %>% sapply(get)
#populations to be renamed and/or merged:

#Filter pop. If S4 genind objects recognizes that the number of rows in each
#list whithin @other, it will filter the correspoding samples in
#@other$wathever_list accordindly. This will not happen if nrows in that given
#list is different to nInd.

for (i in seq_along(gen)){
  if (any(rm_pop[[i]] %in% pop(gen[[i]])))
    gen[[i]] <- gen[[i]][!pop(gen[[i]]) %in% rm_pop[[i]]]
  assertthat::assert_that(nrow(gen[[i]]$other$metadata) == nInd(gen[[i]]))
  }
rm(i)

#assert all of the pop names == metadata$locality names
gen %>% sapply(function(x){
  assertthat::assert_that(
    all(
      x$other$metadata$locality == pop(x)
      ))})

#Change population names from metadata so as to match localities between
# datasets
#Hyla:
# change names
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
# merge populations
#  1
# gen$dart_hyla$other$metadata$locality %<>%
#   gsub(pattern = "Peña de Francia",
#     replacement = "Peña de Francia - Ciudad Rodrigo")
# gen$usat_hyla$other$metadata$locality %<>%
#   gsub(pattern = "Ciudad Rodrigo",
#     replacement = "Peña de Francia - Ciudad Rodrigo")
#   2
gen$dart_hyla$other$metadata$locality %<>%
  gsub(pattern = "Saceruela|Río Esteras",
    replacement = "Sierra de Saceruela")
gen$usat_hyla$other$metadata$locality %<>%
  gsub(pattern = "Fontanosas", replacement = "Sierra de Saceruela")
# 3
gen$usat_hyla$other$metadata$locality %<>%
  gsub(pattern = "Los Baños de Robledillo|Navas de Estena",
    replacement = "Navas de Estena - Baños Robledillo")
gen$dart_hyla$other$metadata$locality %<>%
  gsub(pattern = "Navas de Estena",
    replacement = "Navas de Estena - Baños Robledillo")
#   4
gen$dart_hyla$other$metadata$locality %<>%
  gsub(pattern = "Codesal|Ferreras de Arribas",
    replacement = "Codesal - Ferreras de Arriba")
gen$usat_hyla$other$metadata$locality %<>%
  gsub(pattern = "Codesal|Ferreras de Arriba",
    replacement = "Codesal - Ferreras de Arriba")
#Pelobates:
# change names
gen$dart_pelo$other$metadata$locality %<>%
  gsub(pattern = "Plage Biscarrose", replacement = "Plage Biscarrosse") %>%
  gsub(pattern = "La Porroquia", replacement = "La Parroquia")
# merge population names
# 1
gen$usat_pelo$other$metadata$locality %<>%
  gsub(pattern = "Benalup-Casas Viejas", replacement = "Los Alcornocales")
gen$dart_pelo$other$metadata$locality %<>%
  gsub(pattern = "Tarifa", replacement = "Los Alcornocales")
##END edit pop names
#Assert localities in dartseq and usat are the same
assertthat::assert_that(
  #Pelobates
  identical(
  levels(as.factor(gen$usat_hyla$other$metadata$locality)),
  levels(as.factor(gen$dart_hyla$other$metadata$locality))) &
  #Hyla
  identical(
    levels(as.factor(gen$usat_pelo$other$metadata$locality)),
    levels(as.factor(gen$dart_pelo$other$metadata$locality)))
  )

########## 2nd part ##########

#Keep only shared individuals for Pelobates
# remove individuals from genind obj
# vector with individuals to keep
samples_overlapped_dartusat <-
  indNames(gen$dart_pelo)[(adegenet::indNames(gen$dart_pelo) %in% adegenet::indNames(gen$usat_pelo))]
# new DART genind object
gen[["dart_pelo"]] <-
  gen$dart_pelo[row.names(gen$dart_pelo@tab) %in% samples_overlapped_dartusat]
# new usat genind object
gen[["usat_pelo"]] <-
  gen$usat_pelo[row.names(gen$usat_pelo@tab) %in% samples_overlapped_dartusat]

# check number of Individuals is
assertthat::assert_that(nInd(gen[["usat_pelo"]]) == nInd(gen[["dart_pelo"]]))


#consolite pop slot
gen %<>% lapply(function(x){
  #assert metadata$sample_id are the same to indNames
  assertthat::assert_that(
    all(x$other$metadata$sample_id == indNames(x)))
  pop(x) <- x$other$metadata$locality
  assertthat::assert_that(identical(x$other$metadata$locality,
    as.character(pop(x))))
  x
  })
plot_samples(meta = gen,
  name_plot = "data/intermediate/post_sample_fil_distribution.pdf", res = 5000)
#Save object
saveRDS(gen, "data/intermediate/gen_consolidated_filtered.rds")
