#takes genlight/genind object. Removes samples in metadata non present in
#indNames and reorders metadata$sample_id according to order in indNames
consolidate_metadata <- function(genotypes){
  #ASSERT that
  #1. class is genlight or genind
  assertthat::assert_that(
    class(genotypes) == "genlight" | class(genotypes) == "genind")
  #2. metadata is a sf object
  assertthat::assert_that(class(genotypes@other)[1] == "sf")
  #3. Non duplicated IDs in indNames or metadata
  assertthat::assert_that(
  list(adegenet::indNames(genotypes), genotypes$other$sample_id) %>%
    sapply(function(x) all(duplicated(x)) == F) %>% all())
  #4. All indNames are present in metadata
  assertthat::assert_that(
    all(adegenet::indNames(genotypes) %in% genotypes$other$sample_id))
  #CONSOLIDATE metadata
  #1. remove samples from metadata that are not in indNames
  genotypes$other %<>%
    filter(sample_id %in% indNames(genotypes)) %>%
  #2. reorder rows in metadata according to order in indNames
    arrange(match(sample_id, indNames(genotypes)))
  #assert metadata$sample_id is the same to indNames
  assertthat::assert_that(
    identical(as.character(genotypes$other$sample_id),
              adegenet::indNames(genotypes)))
  genotypes
}
