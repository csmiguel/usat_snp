#create raster file for the Iberian Peninsula
path <- "data/raw/raster_iberian_peninsula.grd"
assertthat::assert_that(!file.exists(path),
    msg = paste(path, "already exists"))

mid_long <- c(-7.5, -2.5, 2.5)
mid_lat <- c(37.5, 42.5)
mid_latlon <- expand.grid(mid_lat, mid_long)

o <- as.list(numeric(nrow(mid_latlon)))
#use SRTM 90 with 90 m resolution at the equator.
#for tiles see:http://srtm.csi.cgiar.org/srtmdata/
for (i in seq_along(o)){
  o[[i]] <- raster::getData("SRTM", download = T, "data/intermediate",
                            lat = mid_latlon[i, 1], lon = mid_latlon[i, 2]) %>%
    raster::aggregate(fact = 3)
}

o$fun <- mean
o$na.rm <- TRUE
mapr <- do.call(raster::mosaic, o)
#saveRDS(mapr, path) would not work because it creates a temporary files
#in the computer which is later removed. Instead it is necessary to writeRaster
raster::writeRaster(mapr, filename = path, format = "raster")
