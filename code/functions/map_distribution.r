#create map with the distribution of the samples and native distribution of the species
#read raster and shape files
p_raster <- "data/raw/raster_iberian_peninsula.rds"
mapr <- readRDS(p_raster)
  #IUCN distribution
p_pelo <- "data/raw/redlist_species_data_44f009a9-b3b2-4608-a818-24c3dd714900"
p_hyla <- "data/raw/redlist_species_data_eff4cf11-63bb-4131-bba3-104c3358697b"
shp_pelo <- rgdal::readOGR(dsn = p_pelo, layer = "data_0") %>%
  sp::spTransform(mapr@crs)
shp_hyla <- rgdal::readOGR(dsn = p_hyla, layer = "data_0") %>%
  sp::spTransform(mapr@crs)

  #dartseq locations
  p_dartmeta <- ""
  dartmeta <- readRDS()
  #microsatellite locations
  p_gen <- ""
  gen <- readRDS()

assertthat::assert_that(as.character(shp_hyla@data$BINOMIAL)[1] == "Hyla arborea")
assertthat::assert_that(as.character(shp_pelo@data$BINOMIAL)[1] == "Pelobates cultripes")

#Plot distribution
plot(mapr, breaks = c(0, 500, 1000, 1500, 2000, 3500),
      col = gray.colors(5, start = 0.9, end = 0.3)) #issue margins
plot(shp_pelo, add = T, border = "black",
  col = grey(0.5, alpha = 0.3))
title(shp_pelo@data$BINOMIAL[1])

#several options for plot:
# plot(h, col = gray.colors(4))
# plot(h, col = gray.colors(4, start = 0.9, end = 0.3))
# plot(h, breaks = c(0, 500, 1000, 1500, 2000, 3500),
#      col = c("lightgreen", "darkgreen", "yellow", "brown", "purple"))
# plot(h, breaks = c(0, 500, 1000, 1500, 2000, 3500),
#      col = gray.colors(5, start = 0.9, end = 0.3),
#      xlim = c(-10, 5))
# spplot(mapr, maxpixels= 500000,
#     col.regions = gray.colors(5, start = 0.9, end = 0.3),
#     at = c(0, 500, 1000, 1500, 2000, 3500))
