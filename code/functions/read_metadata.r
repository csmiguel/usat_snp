#read metadata for dartseq samples
metadata_Dartseq_Hyla <- function(x){
  f <- "data/raw/metadata_hyla_dartseq.xlsx"
  m <- readxl::read_excel(f)
  m_darseq_id <- m[["dartseqID"]]
  assertthat::assert_that(class(x$hyla) == "genlight")
  g_dartseq_id <- indNames(x$hyla) %>% gsub(" ", "", .)
assertthat::assert_that(!is.null(c(m[["dartseqID"]],
      m[["locality"]],
      m[["id"]])))
  #fix sample names
  indNames(x$hyla) <- sapply(g_dartseq_id, function(x){
    which(m[["dartseqID"]] == x) %>% m[["id"]][.]
  })
  #fix pop names
  pop(x$hyla) <- sapply(g_dartseq_id, function(x){
    which(m[["dartseqID"]] == x) %>% m[["locality"]][.]
  })
  metadata_darthyla <<- m %>%
              select(id, locality, latitude, longitude) %>%
              mutate(latitude = as.numeric(latitude)) %>%
              mutate(longitude = as.numeric(longitude)) %>%
              filter(id %in% indNames(x$hyla)) %>%
              rename(sample_id = id) %>%
              #convert to spatial object
              {sp::SpatialPointsDataFrame(data = .[, c(1, 2),
              drop = FALSE], coords = .[, c("latitude", "longitude")])} %>%
              as("sf") %>%
              sf::st_set_crs(4326)
  assertthat::assert_that(identical(
    x = sort(as.character(indNames(x$hyla))),
    y = sort(as.character(metadata_darthyla$sample_id))))
  return(x)
}

metadata_Dartseq_Pelobates <- function(x){
  f <- "data/raw/2019-01-03_meta_Pelobates.csv"
  m <- read.csv(f)
  m %<>% mutate(locality = gsub(pattern = ",.*$", replacement = "", site)) %>%
    rename(full_locality = site) %>%
    rename(dartseqID = vaucher) %>%
    rename(latitude = lat) %>%
    rename(longitude = lon) %>% dplyr::as_tibble() %>%
    {sp::SpatialPointsDataFrame(data = .[, c(1, 6), drop = FALSE],
                                coords = .[, c("latitude", "longitude")])} %>%
    as("sf") %>%
    sf::st_set_crs(4326)

  assertthat::assert_that(class(x$pelo) == "genlight")
  m_darseq_id <- m[["dartseqID"]]
  g_dartseq_id <- indNames(x$pelo)
  assertthat::assert_that(all(g_dartseq_id %in% m_darseq_id), msg =
  "some names in genotypes are missing in metadata")
  #fix pop names
  assertthat::assert_that(!is.null(c(m[["dartseqID"]], m[["locality"]])))
  pop(x$pelo) <- sapply(g_dartseq_id, function(x){
    which(m[["dartseqID"]] == x) %>% m[["locality"]][.]
  }) %>% unlist()
  metadata_dartpelo <<- m %>% #select only samples present in the genlight object
                filter(dartseqID %in% indNames(x$pelo)) %>%
                rename(sample_id = dartseqID)
  return(x)
}
