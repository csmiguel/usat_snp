#this piece of code has been copied from Jeff's. I have modified at likes marked
#with #m
read_Pelobates_cultripes <- function() {
  ## read data
  f <- "data/raw/PcultripesPhylogeography.gen" #m
  o <- adegenet::read.genepop(f, ncode = 3)
  s <-
    read.table("data/raw/table-s1-data.csv",
               header = TRUE, sep = ",", stringsAsFactors = FALSE) %>%
    tibble::as_tibble() %>%
    dplyr::filter(N > 0) %>%
    plyr::ddply("ID", function(x) {
      # iterate over codes and interpolate codes specified as ranges
      plyr::ldply(
        stringr::str_squish(strsplit(x$Sample.Codes, ",")[[1]]),
        function(y) {
          # extract codes
          if (grepl("-", y, fixed = TRUE)) {
            ## find start and end codes
            start <- strsplit(y, "-", fixed = TRUE)[[1]]
            end <- start[[2]]
            start <- start[[1]]
            ## find character prefix for code
            prefix <- gsub("[[:digit:]]", "", start)
            ## interpolate numbers
            start_number <- as.numeric(gsub(prefix, "", start, fixed = TRUE))
            end_number <- as.numeric(gsub(prefix, "", end, fixed = TRUE))
            assertthat::assert_that(is.finite(start_number),
                                    is.finite(end_number))
            ## create codes
            codes <- paste0(prefix, seq(start_number, end_number))
          } else {
            codes <- y
          }
          # return table with sample codes per locality
          data.frame(ID = x$ID, Locality = x$Locality, Sample.Code = codes,
                     Longitude = x$Longitude, Latitude = x$Latitude,
                     stringsAsFactors = FALSE)
        })
    }) %>%
    tibble::as_tibble()
  ## add in population data,
  ### since we only have sample ids for individuals with mtDNA data,
  ### we don't have complete information on which samples belong to which
  ### populations, but if we known that sample A is in the same population as
  ### as sample B, and we known that sample B is in population A then we know
  ### that sample A is in population A
  s2 <- tibble::tibble(sample_id = rownames(o@tab),
                       population_id = as.character(o@pop)) %>%
       dplyr::left_join(s, by = c("sample_id" = "Sample.Code")) %>%
       plyr::ddply("population_id", function(x) {
         if (!all(is.na(x$ID))) {
           x$ID <- na.omit(x$ID)[[1]]
           x$Locality <- na.omit(x$Locality)[[1]]
           x$Longitude <- na.omit(x$Longitude)[[1]]
           x$Latitude <- na.omit(x$Latitude)[[1]]
         }
         x
       }) %>%
       tibble::as_tibble() %>%
       dplyr::right_join(tibble::tibble(sample_id = rownames(o@tab)),
                         by = "sample_id")
  ### fortunately it appears that only one population did not have any mtDNA
  ### samples so we can manually define the  popuation ids for this population
  assertthat::assert_that(dplyr::n_distinct(s2$population_id[is.na(s2$ID)]) ==
                          1)
  missing_population <- setdiff(unique(s$ID), unique(s2$ID[!is.na(s2$ID)]))
  assertthat::assert_that(assertthat::is.number(missing_population))
  na_pos <- which(is.na(s2$ID))
  repl_pos <- which(s$ID == missing_population)[1]
  s2$ID[na_pos] <- missing_population
  s2$Locality[na_pos] <- s$Locality[repl_pos]
  s2$Longitude[na_pos] <- s$Longitude[repl_pos]
  s2$Latitude[na_pos] <- s$Latitude[repl_pos]
  ## assert that no NA values
  assertthat::assert_that(nrow(s2) == nrow(na.omit(s2)))
  ## manually recode populations
  assertthat::assert_that(identical(s2$sample_id, rownames(o@tab)))
  o@pop <- factor(s2$Locality)
  # add in spatial data
  o@other <-
    s2 %>%
    dplyr::select(sample_id, Locality, Longitude, Latitude) %>%
    dplyr::rename(population_id = Locality, x = Longitude, y = Latitude) %>%
    {sp::SpatialPointsDataFrame(data = .[, c(1, 2), drop = FALSE],
                                coords = .[, c("x", "y")])} %>%
    as("sf") %>%
    sf::st_set_crs(4326)
  ## return result
  o
}

read_Hyla_molleri <- function() {
  ## read data
  g <- read.table("data/raw/appendix-s2.txt",
                  header = FALSE, sep = ";", stringsAsFactors = FALSE,
                  skip = 1, comment.char = "") %>%
       tibble::as_tibble() %>%
       dplyr::rename(sample_id = V1, population_id = V2)
       g[["sample_id"]] <- gsub(pattern = "^.*#", #m
                      replacement = "", g[["sample_id"]]) #m

  s <- read.table("data/raw/table-1-data.csv",
                  header = TRUE, sep = ",", stringsAsFactors = FALSE) %>%
       tibble::as_tibble()
  ## create genind object
  o <- as.matrix(g[, c(-1, -2)])
  mode(o) <- "character"
  o[o == "0"] <- "000"
  o[nchar(o) == 2] <- paste0("0", o[nchar(o) == 2])
  o2 <- matrix("", ncol = ncol(o) / 2, nrow = nrow(o),
               dimnames = list(g$sample_id, colnames(o)[seq(1, ncol(o), 2)]))
  cols <- seq(1, ncol(o), 2)
  for (i in seq_along(cols)){
    o2[, i] <- paste0(o[, cols[i]], "/", o[, cols[i] + 1])
  }
   o2 <- adegenet::df2genind(o2, ploidy = 2, sep = "/", type = "codom",
                            NA.char = "000")
  ## add in population data
  g2 <- dplyr::left_join(g, s, by = c("population_id" = "id"))
  o2@pop <- factor(g2$locality)
  ## add in sample location data
  s2 <- g2[, c("sample_id", "locality", "long", "lat")] %>%
        setNames(c("sample_id", "population_id", "x", "y")) %>%
        dplyr::mutate(x = gsub("−", "-", x, fixed = TRUE)) %>%
        dplyr::mutate(x = as.numeric(x), y = as.numeric(y)) %>%
        {sp::SpatialPointsDataFrame(data = .[, c(1, 2), drop = FALSE],
                                    coords = .[, c("x", "y")])} %>%
        as("sf") %>%
        sf::st_set_crs(4326)
  o2@other <- s2
  ## return data
  o2
}
