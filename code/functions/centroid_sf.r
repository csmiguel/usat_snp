#calculate centroid from a sf object given an index (default = "locality")
centroid_sf <- function(sf_object = sf_object, index = "locality",
                        latlon = c("longitude", "latitude")){
  sfc_as_cols(sf_object) %>%
  plyr::ddply(index, function(x){
  x <- dplyr::select(x, longitude, latitude)
  x <- sp::SpatialPoints(x, proj4string = CRS("+init=epsg:4326"))
  x <- rgeos::gCentroid(x)
  x <- as.data.frame(x)
  names(x) <- latlon
  x
  }) %>%
  {sp::SpatialPointsDataFrame(data = .[, 1, drop = FALSE],
        coords = .[, latlon])} %>%
  as("sf") %>%
  sf::st_set_crs(4326)
}


sfc_as_cols <- function(x, names = c("longitude", "latitude")) {
  stopifnot(inherits(x, "sf") && inherits(sf::st_geometry(x), "sfc_POINT"))
  ret <- sf::st_coordinates(x)
  ret <- tibble::as_tibble(ret)
  stopifnot(length(names) == ncol(ret))
  x <- x[, !names(x) %in% names]
  ret <- setNames(ret, names)
  x <- dplyr::bind_cols(x, ret)
  st_geometry(x) <- NULL
  x
}
