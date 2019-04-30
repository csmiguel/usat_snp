#read metadata hyla
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
  metadata_hyla <<- m
  return(x)
}

metadata_Dartseq_Pelobates <- function(x){
  f <- "data/raw/2019-01-03_meta_Pelobates.csv"
  m <- read.csv(f)
  m %<>% mutate(locality = gsub(pattern = ",.*$", replacement = "", site)) %>%
    rename(full_locality = site) %>%
    rename(dartseqID = vaucher) %>%
    rename(latitude = lat) %>%
    rename(longitude = lon)

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
  metadata_pelo <<- m
  return(x)
}
