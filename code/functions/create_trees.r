#create trees
# gen is a genlight or genind of object
create_trees <- function(gen, type = c("usat", "snp")){
  source("code/functions/manhattan_dist.r")
assertthat::assert_that(class(gen) == "genind", msg = "not genind")
  if (type == "usat"){
    assertthat::assert_that(adegenet::nLoc(gen) < 100, msg = "too many loci")
    gen.d <- genind_manhattan(gen)
  }
  else if (type == "snp"){
    assertthat::assert_that(adegenet::nLoc(gen) > 1000, msg = "too few loci")
    if (class(gen) == "genind")
    gen.d <- stats::dist(gen, method = "manhattan") / 2
  } else (stop("define genotype format"))
  z <- ape::njs(gen.d) #issue2: is nj equal njs
  z
 }
