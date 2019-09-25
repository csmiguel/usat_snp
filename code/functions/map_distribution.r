plot_samples <- function(meta = metadata, name_plot, res = res,
  mode = c(1, 2, 3)){
#meta is a metadata list where the name of each element has names of marker
#and species. Each element is a sf object with #sample_id#locality#geometry
#create map with the distribution of the samples and native
#distribution of the species
#res, resolution of raster. For fast plotting use 5.10^3. For high-res 5.10^4
#mode is 1 for test mode and 2 for pretty plotting; 3 is for plotting #id
library(dplyr)
library(sf)
library(raster)
assertthat::assert_that(seq_along(meta) %>% sapply(function(x){
  names(meta[[1]]@other$metadata) == c("sample_id", "locality",  "geometry")
  }) %>% all())
#read raster and shape files
  #raster
p_raster <- "data/raw/raster_iberian_peninsula.grd"
mapr <- raster::raster(p_raster)
  #IUCN distribution
p_pelo <- "data/raw/redlist_species_data_44f009a9-b3b2-4608-a818-24c3dd714900"
p_hyla <- "data/raw/redlist_species_data_eff4cf11-63bb-4131-bba3-104c3358697b"
shp_pelo <- rgdal::readOGR(dsn = p_pelo, layer = "data_0") %>%
  sp::spTransform(mapr@crs) %>% raster::crop(extent(mapr))
shp_hyla <- rgdal::readOGR(dsn = p_hyla, layer = "data_0") %>%
  sp::spTransform(mapr@crs) %>% raster::crop(extent(mapr))
    #check I loaded the right layers
assertthat::assert_that(as.character(
  shp_hyla@data$BINOMIAL)[1] == "Hyla arborea")
assertthat::assert_that(as.character(
  shp_pelo@data$BINOMIAL)[1] == "Pelobates cultripes")

#load functions
source("code/functions/centroid_sf.r")

#Plot distribution
outfile <- name_plot
ratio <- 2
size <- 5
x <- c("hyla", "pelo")
y <- c("Hyla molleri", "Pelobates cultripes")
micro_col <- "blue"
snp_col <- "red"
pch_snp <- 20
pch_usat <- 2
cex_pt <- c(1, 1.4)

pdf(outfile, width = size * ratio, height = size)
par(mfrow = c(1, 2))
for (i in seq_along(x)){
  #sf centroids usats
  centroids_usat <- centroid_sf(
    meta[[grep(paste0("usat.", x[i]), names(gen))]]@other$metadata)
  #sf centroids dart
  centroids_dart <- centroid_sf(
    meta[[grep(paste0("dart.", x[i]), names(gen))]]@other$metadata)
  #plot
  plot(mapr, maxpixels = res, #raster
       breaks = c(0, 500, 1000, 1500, 2000, 3500),
       col = gray.colors(5, start = 0.9, end = 0.3),
       axes = T, box = FALSE, ylab = "Latitude", xlab = "Longitude")
  plot(get(ls()[grepl(paste0("shp_", x[i]), ls())]),  #IUCN distribution
       add = T, border = "black", col = grey(0.5, alpha = 0.3))
  plot(sf::st_geometry(centroids_dart), add = T, pch = pch_snp, cex = cex_pt[2],
       col = snp_col) #snps samples
  plot(sf::st_geometry(centroids_usat), add = T, pch = pch_usat, cex = cex_pt[1],
       col = micro_col) #usats samples
  title(y[i], font.main = 3, cex.main = 0.8)
  if (i == 1){
    legend(y = extent(mapr)@ymax, x = extent(mapr)@xmin,
      legend = c("Microsatellites", "SNPs"), pt.cex = cex_pt, cex = 0.7,
      pch = c(pch_usat, pch_snp), bty = "n", col = c(micro_col, snp_col))
  }
  off <- 0.5
  pop_usat <- centroids_usat %>% sfc_as_cols %>%
    mutate(longitude = as.numeric(longitude)) %>%
    mutate(latitude = as.numeric(latitude))
  pop_dart <- centroids_dart %>% sfc_as_cols %>%
    mutate(longitude = as.numeric(longitude)) %>%
    mutate(latitude = as.numeric(latitude))
    if (mode == 1){
  text(pop_usat$longitude, y = pop_usat$latitude, labels = pop_usat$locality,
    cex = .2, pos = 1, offset = .2, col = micro_col, font = 2)
  text(pop_dart$longitude, y = pop_dart$latitude, labels = pop_dart$locality,
    cex = .2, pos = 3, offset = .2, col = snp_col, font = 2)
  }
  if (mode == 2){
    text(pop_usat$longitude, y = pop_usat$latitude, labels = pop_usat$locality,
      cex = .3, pos = 1, offset = .3, col = "black", font = 2)
  }
  if (mode == 3){
    t1 <- readRDS("data/intermediate/table1.rds")
    text(pop_usat$longitude, y = pop_usat$latitude,
      labels = paste(t1$ID[match(pop_usat$locality, t1$locality)],
                     pop_usat$locality, sep = "-"),
      cex = .2, pos = 1, offset = .3, col = "black", font = 2)
  }
}
dev.off()
}
