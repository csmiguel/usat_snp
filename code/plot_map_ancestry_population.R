###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: Plot map with ancestries per population
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

library(dplyr)
library(magrittr)
library(sf)
library(raster)
library(adegenet)

#vector to reorder datasets
dor <- c("dart_hyla", "usat_hyla", "dart_pelo", "usat_pelo")

#genotypes
gen <- readRDS("data/intermediate/gen_consolidated_filtered.rds")
gen <- gen[match(dor, names(gen))]

#ancestries for major clusters
qclumpak <- readRDS("data/intermediate/clumpak_major.rds")
# remove subset of loci
qclumpak %<>% .[grep("^[a-z]", names(qclumpak))] %>% {.[match(dor, names(.))]}

#load plotting function
source("code/functions/map_distribution.r")
#load functions
source("code/functions/centroid_sf.r")
#colors for ancestries
source("code/parameters/structure_colors.r")

#assertions

assertthat::assert_that(seq_along(gen) %>% sapply(function(x){
  names(gen[[1]]@other$metadata) == c("sample_id", "locality",  "geometry")
}) %>% all())

#1. Read raster and shape files
#raster
p_raster <- "data/raw/raster_iberian_peninsula.grd"
mapr <- raster::raster(p_raster)
#IUCN distribution
p_pelo <- "data/raw/redlist_species_data_44f009a9-b3b2-4608-a818-24c3dd714900"
p_hyla <- "data/raw/redlist_species_data_eff4cf11-63bb-4131-bba3-104c3358697b"
shp_pelo <- rgdal::readOGR(dsn = p_pelo, layer = "data_0") %>%
  sp::spTransform(mapr@crs)
shp_hyla <- rgdal::readOGR(dsn = p_hyla, layer = "data_0") %>%
  sp::spTransform(mapr@crs)
#check I loaded the right layers
assertthat::assert_that(as.character(
  shp_hyla@data$BINOMIAL)[1] == "Hyla arborea")
assertthat::assert_that(as.character(
  shp_pelo@data$BINOMIAL)[1] == "Pelobates cultripes")


#compute centroids for each metadata
centroids <-
  seq_along(gen) %>%
  lapply(function(x){centroid_sf(gen[[x]]@other$metadata)})
names(centroids) <- names(gen)

#compute mean ancestries per location
# it returns a list with the same structure as qclumpak but with mean ancestries per location
assertthat::assert_that(all(names(gen) == names(qclumpak)))
mean_anc <-
  seq_along(gen) %>%
  lapply(function(x){
    qlist <- qclumpak[[x]]
    meta <- gen[[x]]@other$metadata
    h <-
      2:length(qlist) %>%
      lapply(function(y){
        split(qlist[[y]], meta$locality) %>%
        lapply(function(s) apply(s, 2, mean)) %>%
          do.call(what = rbind)
        })
    names(h) <- names(qclumpak[[1]])[2:length(qlist)]
    h
    })
names(mean_anc) <- names(gen)

#plot
reso = 50000
pdf(file = "data/final/maps_ancestry.pdf",
    height = 11, width = 8)
par(mfrow = c(7, 4), mar = c(0.5, 0, 0.5, 0), lwd = 0.2)
for (k in seq_along(mean_anc[[1]])){
  for (dataset in seq_along(gen)){
    plot(mapr, maxpixels = reso,
         breaks = c(0, 500, 1000, 1500, 2000, 3500),
         col = gray.colors(5, start = 0.9, end = 0.4),
         axes = F, box = FALSE, legend = FALSE)
    #get IUCN shape file
    shpf <-
    get(ls()[match(names(gen)[dataset] %>%
                     sub(pattern = "^.*_", replacement = "shp_"), ls())])
    plot(shpf, add = T, border = "black", col = grey(0.7, alpha = 0.3))
    #overplot ancestries as pieplots
    coo <- centroids[[names(gen)[dataset]]] %>% sfc_as_cols()
    anc <- mean_anc[[names(gen)[dataset]]][[k]]
    mapplots::draw.pie(x = coo$longitude,
                       y = coo$latitude,
                       z = anc,
                       radius = 0.25, lty = 1, col = str_col[1:ncol(anc)])
  }
}
dev.off()
