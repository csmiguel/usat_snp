#bootstrap trees
   #nboot: number of bootstrap
   #type: "usat" "snp"; usat is a genind object with pop slot filled with pops
   #     "snp" is a genlight object
bootstrap_nj <- function(gen, nboot = nboot, type = c("usat", "snp")){
  #assertions and conversions
  source("code/functions/manhattan_dist.r")
  assertthat::assert_that(class(gen) %in% c("genind", "genlight"),
    msg = "not genind or genlight")
  if (type == "usat"){
    if(adegenet::nLoc(gen) < 100) warning("too many loci")
    gen_matrix <- adegenet::genind2df(gen) %>% .[, -1] #convert genind 2 df
  }
  if (type == "snp"){
    if(adegenet::nLoc(gen) > 100) warning("too few loci")
    if (class(gen) == "genind") gen <- dartR::gi2gl(gen)
    gen_matrix <- as.matrix(gen)
  }
  n_col <- adegenet::nLoc(gen)
  zz <- list()
  #begin bootstrapping
  for (i in 1:nboot){
    #calculate distance matrices
    if (type == "usat"){
      gen.d <-
        gen_matrix[, sample(1:n_col, size = n_col, replace = T)] %>%
        setNames(1:n_col) %>%
        adegenet::df2genind(ncode = 3, type = "codom") %>%
        genind_manhattan()
  }
  else if (type == "snp"){
    gen.d <-
      gen_matrix[, sample(1:n_col, size = n_col, replace = T)] %>%
      stats::dist(method = "manhattan")
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
get_bt <- function(list_tr, external_ref_tree = NULL){
  #list_tr is a list of multiphylo objects
  h <- names(list_tr)
  if (!is.null(external_ref_tree)){
    assertthat::assert_that(
      ape::Ntip(external_ref_tree) == ape::Ntip(list_tr[[1]][1]))
    list_tr <-
      seq_along(list_tr) %>%
      lapply(function(x){
        c(external_ref_tree, list_tr[[x]])
      })
    }
  bt <-
    seq_along(list_tr) %>%
    sapply(function(x){
      temp <- list()
      #1. for rooted trees
      if (all(sapply(list_tr, ape::is.rooted))){
        pp <- ape::prop.part(list_tr[[x]][-1])
        ans <- ape::prop.clades(list_tr[[x]][[1]], part = pp, rooted = T)
      }
      #2. for unrooted trees
      if (all(!sapply(list_tr, ape::is.rooted))){
        phy <- reorder(list_tr[[x]][[1]], "postorder")
        ints <- phy$edge[, 2] > ape::Ntip(phy)
        ans <- ape::countBipartitions(phy, list_tr[[x]][-1])
        ans <- c(nboot, ans[order(phy$edge[ints, 2])])
      }
      ans[is.na(ans)] <- 0
      temp[[1]] <- ans
      temp
    })
    names(bt) <- h
    bt
  }
