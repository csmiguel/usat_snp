#bootstrap trees
   #nboot: number of bootstrap
   #type: "usat" "snp"; usat is a genind object with pop slot filled with pops
   #     "snp" is a genlight object
create_trees <- function(gen, type = c("usat", "snp")){
  if (type == "usat"){
    assertthat::assert_that(class(gen) == "genind", msg = paste0("not genind"))
    assertthat::assert_that(nLoc(gen) < 100, msg = "too many loci")
    gen.d <- hierfstat::genind2hierfstat(gen) %>%
    hierfstat::genet.dist(method = "Da")
  }
  else if (type == "snp"){
    assertthat::assert_that(class(gen) == "genlight", msg = "not genlight")
    assertthat::assert_that(ncol(gen) > 1000, msg = "too few loci")
    gen.d <- as.matrix(gen) %>% dist()
  } else (stop("define genotype format"))
  z <- ape::njs(gen.d) #issue2: is nj equal njs
  z
 }
