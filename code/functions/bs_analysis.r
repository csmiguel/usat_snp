tree_metrics <- function(tr, bs, external_ref_tree = NULL){
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

  if (!is.null(external_ref_tree)){
    tr1 <- external_ref_tree
    } else {
      tr1 <- tr[1][[1]]
      }#reference tree
  bs <- bs / nboot #vector with bs values normalized to 1
  #Normalize length of the tree
  #distance between tips only
  disttips <- ape::cophenetic.phylo(tr1) %>% sort()
  #distance between all tips and nodes, but selecting only the quandrant of the
  # matrix with distances between tips
  distnodes <- ape::dist.nodes(tr1)[1:Ntip(tr1), 1:Ntip(tr1)] %>% sort()
  #asssert that the order of the columns is tips + nodes
  assertthat::assert_that(all(disttips == distnodes),
              msg = paste0("check order of tips and nodes in trees"))
  #meaning that the bottom-right quadrant of the matrix corresponds to
  #distances between internal nodes
  maxd <- ape::dist.nodes(tr1)[(ape::Ntip(tr1) + 1):(tr1$Nnode + ape::Ntip(tr1)),
                             (ape::Ntip(tr1) + 1):(tr1$Nnode + ape::Ntip(tr1))] %>% max()
  # normalize
  tr1$edge.length <- tr1$edge.length / maxd
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
