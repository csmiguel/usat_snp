#read metadata for dartseq samples
metadata_Dartseq_Hyla <- function(x){
  f <- "data/raw/metadata_hyla_dartseq.xlsx"
  m <- readxl::read_excel(f)
  m_darseq_id <- m[["dartseqID"]]
  assertthat::assert_that(class(x$hyla) == "genlight")
  g_dartseq_id <- adegenet::indNames(x$hyla) %>% gsub(" ", "", .)
assertthat::assert_that(!is.null(c(m[["dartseqID"]],
      m[["locality"]],
      m[["id"]])))
  #fix sample names
  adegenet::indNames(x$hyla) <- sapply(g_dartseq_id, function(x){
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
              filter(id %in% adegenet::indNames(x$hyla)) %>%
              rename(sample_id = id) %>%
              #convert to spatial object
              {sp::SpatialPointsDataFrame(data = .[, c(1, 2),
              drop = FALSE], coords = .[, c("longitude", "latitude")])} %>%
              as("sf") %>%
              sf::st_set_crs(4326)
  assertthat::assert_that(identical(
    x = sort(as.character(adegenet::indNames(x$hyla))),
    y = sort(as.character(metadata_darthyla$sample_id))))
  return(x)
}

metadata_Dartseq_Pelobates <- function(x){
  f <- "data/raw/meta_pelobates.csv"
  m <- read.delim(f, sep = ";")
  m %<>%
    mutate(locality = gsub(pattern = ",.*$", replacement = "", site)) %>%
    mutate(collector = gsub(pattern = " ", replacement = "", collector)) %>%
    rename(full_locality = site) %>%
    rename(dartseqID = vaucher) %>%
    rename(latitude = lat) %>%
    rename(longitude = lon) %>% dplyr::as_tibble() %>%
    {sp::SpatialPointsDataFrame(data = .[, c(1, 6, 7), drop = FALSE],
                                coords = .[, c("longitude", "latitude")])} %>%
    as("sf") %>%
    sf::st_set_crs(4326)

  assertthat::assert_that(class(x$pelo) == "genlight")
  m_darseq_id <- m[["dartseqID"]]
  g_dartseq_id <- adegenet::indNames(x$pelo)
  assertthat::assert_that(all(g_dartseq_id %in% m_darseq_id), msg =
  "some names in genotypes are missing in metadata")
  #fix pop names
  assertthat::assert_that(!is.null(c(m[["collector"]], m[["locality"]])))

  pop(x$pelo) <- sapply(g_dartseq_id, function(x){
    which(m[["dartseqID"]] == x) %>% m[["locality"]][.]
  }) %>% unlist()

  indNames(x$pelo) <- sapply(g_dartseq_id, function(x){
    which(m[["dartseqID"]] == x) %>% m[["collector"]][.]
  }) %>% unlist()

  attr(x = pop(x$pelo), "names") <- indNames(x$pelo)
  m %<>% dplyr::select(-dartseqID)
  metadata_dartpelo <<- m %>% #select only samples present in the genlight object
                filter(collector %in% adegenet::indNames(x$pelo)) %>%
                rename(sample_id = collector)
  return(x)
}
