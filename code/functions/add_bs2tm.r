#Add bootstrap values from ape prop.clades to treeman object
add_bs2tm <- function(tr, bs, tm){
  #tr is a phylo tree used to get the support values. it should be the reference
  # rooted tree over which prop.claes are calculated.
  #bs is a vector with bootstrap support values calculated with function
  # ape::prop.par
  #tm is a treeman object obtained from as(tr, "TreeMan")
  #-> the function finds matches between the clades in the treeman tree and the partitions in a phylo tree based on tip names. It then adds a slot in the treeman object with bootstrap support values.
  assertthat::assert_that(length(bs) == tr$Nnode)
  assertthat::assert_that(all(tm@tips %in% tr$tip.label))
  assertthat::assert_that(ape::is.rooted(tr))
  # 1. list of names in partitions of phylo tree
  # get partitions from phylo tree
  tr_p <- ape::prop.part(tr)
  # replace tip numbers by tip names in the list
  tr_pl <-
    seq_along(tr_p) %>% # for each partition
    #replace tip No by tip name in the list
    lapply(function(x){
      plyr::mapvalues(tr_p[[x]],
                      from = 1:ape::Ntip(tr),
                      to = attr(tr_p, "labels"),
                      warn_missing = F)
      })

  #function: for each clade in treeman tree it finds its equivalent partition
  # in the partition (prop.part) object
  bs_a <- function(x){
    #tip names for node x
    tips_ndxtm <- treeman::getSubtree(tm, id = names(tm@ndlst)[x])@tips
    #position of node x in partition list
    pos_ndx <-
      seq_along(tr_pl) %>%
      sapply(function(y){
        identical(sort(tips_ndxtm), sort(tr_pl[[y]]))
      }
      ) %>% which(. == T)
  }

  # 2. vector with position of @ndlst nodes in phylo partitions.
  pos <-
    seq_along(tm@ndlst) %>% #list of nodes in treeman tree
      sapply(function(x){
        #if it is not a node
        if(!names(tm@ndlst)[x] %in% tm@nds)
          a <- NA
        #if it is a node
        if(names(tm@ndlst)[x] %in% tm@nds){
          a <- bs_a(x)
        }
        a
      })
  #3. add new slot with bootstrap support values to treeman
  # node ids from @ndlst
  ids <- names(tm@ndlst) %>% grep(pattern = "n", value = T)
  # bootstrap support values
  vals <- bs[pos[!is.na(pos)]]
  # create new slot with bs values
  treeman::setNdsOther(tree = tm, ids = ids, vals = vals, slt_nm = "bs") %>%
      updateSlts
}
