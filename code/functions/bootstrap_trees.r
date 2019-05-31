#bootstrap trees
   #nboot: number of bootstrap
   #type: "usat" "snp"; usat is a genind object with pop slot filled with pops
   #     "snp" is a genlight object
bootstrap_nj <- function(gen, nboot = nboot, type = c("usat", "snp")){
  #assertions and conversions
  assertthat::assert_that(class(gen) %in% c("genind", "genlight"),
    msg = "not genind nor genlight")
  if (type == "usat"){
    assertthat::assert_that(adegenet::nLoc(gen) < 100, msg = "too many loci")
    genh <- hierfstat::genind2hierfstat(gen)
    ind_n <- adegenet::indNames(gen)
    gen_matrix <- genh[, -1] #genotypes without sample names
  }
  if (type == "snp"){
    assertthat::assert_that(adegenet::nLoc(gen) > 100, msg = "too few loci")
    if (class(gen) == "genind") gen <- dartR::gi2gl(gen)
    gen_matrix <- as.matrix(gen)
  }
  n_col <- ncol(gen_matrix)
  zz <- list()
  #begin bootstrapping
  for (i in 1:nboot){
    #calculate distance matrices
    if (type == "usat"){
      gen.d <- gen_matrix[, sample(1:n_col, size = n_col, replace = T)] %>%
      {hierfstat::genet.dist(data.frame("pop" = ind_n, .), method = "Da")}
      attr(gen.d, "Labels") <- ind_n #add sample names
  }
  else if (type == "snp"){
    gen.d <- gen_matrix[, sample(1:n_col, size = n_col, replace = T)] %>%
    dist()
  } else (stop("define genotype format"))
  #neighbor joining trees
  zz[[i]] <- ape::njs(gen.d)
 }
 class(zz) <- "multiPhylo"
 zz
}

#root trees by midpoint rooting
midpoint_tr <- function(list_mp){
  #list_mp is a list of multiphylo objects
  zz <- list()
  for (i in seq_along(list_mp)){
    h <- lapply(list_mp[[i]], phytools::midpoint.root)
    class(h) <- "multiPhylo"
    zz[[i]] <- h
  }
  names(zz) <- names(list_mp)
  zz
}

#get bootstrap value for nodes for rooted or unrooted trees
get_bt <- function(list_tr){
  #list_tr is a list of multiphylo objects
  bt <-
    seq_along(list_tr) %>%
    sapply(function(x){
      temp <- list()
      #1. for rooted trees
      if (all(sapply(list_tr, ape::is.rooted))){
        pp <- ape::prop.part(list_tr[[x]][-1])
        ans <- ape::prop.clades(list_tr[[x]][[1]], part = pp, rooted = T)
      }
      #1. for unrooted trees
      if (all(!sapply(list_tr, ape::is.rooted))){
        phy <- reorder(list_tr[[x]][[1]], "postorder")
        ints <- phy$edge[, 2] > ape::Ntip(phy)
        ans <- ape::countBipartitions(phy, list_tr[[x]][-1])
        ans <- c(nboot, ans[order(phy$edge[ints, 2])])
      }
      temp[[1]] <- ans
      temp
    })
    names(bt) <- names(list_tr)
    bt
  }
