#create map with the distribution of the samples and native distribution of the species
library(dplyr)
library(sf)
library(raster)
#read raster and shape files
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

  #sample locations
  p_meta <- "data/intermediate/metadata_all.rds"
  meta <- readRDS(p_meta)
#load functions
source("code/functions/centroid_sf.r")

#Plot distribution
outfile <- paste0("data/final/map_samples.pdf")
  ratio <- 2
  size <- 5
x <- c("hyla", "pelo")
y <- c("Hyla molleri", "Pelobates cultripes")
micro_col <- "blue"
snp_col <- "red"
pch_snp <- 20
pch_usat <- 2

pdf(outfile, width = size * ratio, height = size)
par(mfrow = c(1, 2))
for (i in seq_along(x)){
  plot(mapr, maxpixels = 50000,
       breaks = c(0, 500, 1000, 1500, 2000, 3500),
       col = gray.colors(5, start = 0.9, end = 0.3),
       axes = T, box = FALSE)
  plot(get(ls()[grepl(paste0("shp_", x[i]), ls())]),  #IUCN
       add = T, border = "black", col = grey(0.5, alpha = 0.3))
  plot(sf::st_geometry(centroid_sf(meta[grepl(x[i], names(meta))][[1]])),
       add = T, pch = pch_snp, cex = 1.4, col = snp_col) #snps
  plot(sf::st_geometry(centroid_sf(meta[grepl(x[i], names(meta))][[2]])),
       add = T, pch = pch_usat, cex = 1, col = micro_col) #usats
  title(y[i], font = 3)
  if (i == 1){
    legend("topleft", legend = c("Microsatellites", "SNPs"),
           pch = c(pch_usat, pch_snp), bty = "n", col = c(micro_col, snp_col))
  }
}
dev.off()
