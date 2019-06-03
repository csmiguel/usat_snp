tree_metrics <- function(tr, bs){
  #tr is a multiphylo object
  #bs is a vector with bootstrap supports
  #returns a data frame with
  # dnp: distance from a node to its parent node
  # dnr: distance from a node to the root
  # bs_tm: bootstrap support
  library(ggplot2)
  require(gridExtra)
  library(treeman)
  source("code/functions/add_bs2tm.r")
  #
  tr1 <- tr[1][[1]]#reference tree
  bs <- bs / nboot #vector with bs values normalized to 1
  #Normalize length of the tree
  # distance between all nodes (also tips) in the tree
  maxd <- ape::cophenetic.phylo(tr1) %>% max()
  # normalize
  tr1$edge.length <- tr1$edge.length / maxd
  #assert the max distance is 1 => normalization is correct
  assertthat::assert_that(round(max(ape::cophenetic.phylo(tr1)), 3) == 1)
  #Convert to treeman object
  tm <- as(tr1, "TreeMan")
  # add bootstrap support to treeman
  tm <- add_bs2tm(tr = tr1, bs = bs, tm = tm)
  #metrics to summarize the tree:
  # 1. distance from a node to its parent node
  dnp <-
    tm@ndlst[names(tm@ndlst) %in% tm@nds] %>%
      sapply(function(x){
        x$spn
      })
  #   2. distance from nodes to root
  dnr <-
    seq_along(tm@nds) %>%
      sapply(function(x){
        treeman::getNdPrdst(tm, tm@nds[x])
      })
  names(dnr) <- tm@nds
  #reorder names of dnr to match values of dnp
  dnr <- dnr[match(names(dnp), names(dnr))]
  # 3. extract sorted support values
  bs_tm <- tm["bs"][!is.na(tm["bs"])]
  #bind to data frame
  rbind(dnp, dnr, bs_tm) %>%
   t %>%
   tibble::as_tibble() %>%
   dplyr::filter(dnr > 0) #remove root
}
