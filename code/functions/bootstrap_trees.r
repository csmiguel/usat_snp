#bootstrap trees
   #nboot: number of bootstrap
   #type: "usat" "snp"; usat is a genind object with pop slot filled with pops
   #     "snp" is a genlight object
bootstrap_nj <- function(gen, nboot = 100, type = c("usat", "snp")){
  zz <- list()
  for (i in 1:nboot){
    if (type == "usat"){
      assertthat::assert_that(class(gen) == "genind", msg = paste0("not genind"))
      assertthat::assert_that(nLoc(gen) < 100, msg = "too many loci")
      genh <- hierfstat::genind2hierfstat(gen)
      gen_matrix <- genh[, -1] #genotypes without sample names
      gen.d <- gen_matrix[, sample(1:ncol(gen_matrix),
        size = ncol(gen_matrix), replace = T)] %>%
      {hierfstat::genet.dist(data.frame("pop" = genh[, 1], .), method = "Da")}
  }
  else if (type == "snp"){
    assertthat::assert_that(class(gen) == "genlight", msg = "not genlight")
    assertthat::assert_that(ncol(gen) > 100, msg = "too few loci")
    gen_matrix <- as.matrix(gen)
    gen.d <- gen_matrix[, sample(1:ncol(gen_matrix),
      size = ncol(gen_matrix), replace = T)] %>% dist()
  } else (stop("define genotype format"))
  zz[[i]] <- ape::njs(gen.d) #issue2: is nj equal njs
 }
 class(zz) <- "multiPhylo"
 zz
}
