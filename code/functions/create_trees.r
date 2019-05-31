#create trees
# gen is a genlight or genind of object
create_trees <- function(gen, type = c("usat", "snp")){
assertthat::assert_that(class(gen) %in% c("genlight", "genind"),
  msg = "not genlight or genind")
  if (type == "usat"){
    assertthat::assert_that(adegenet::nLoc(gen) < 100, msg = "too many loci")
    ind_n <- adegenet::indNames(gen)
    gen.d <- hierfstat::genind2hierfstat(gen) %>%
      dplyr::mutate(pop = ind_n) %>%
      hierfstat::genet.dist(method = "Da")
    attr(gen.d, "Labels") <- ind_n
  }
  else if (type == "snp"){
    assertthat::assert_that(adegenet::nLoc(gen) > 1000, msg = "too few loci")
    if (class(gen) == "genind")
    gen <- dartR::gi2gl(gen)
    gen.d <- as.matrix(gen) %>% dist()
  } else (stop("define genotype format"))
  z <- ape::njs(gen.d) #issue2: is nj equal njs
  z
 }
