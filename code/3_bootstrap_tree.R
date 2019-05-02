###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: produce phylogenetic trees for microsatellite and snp data
#DESCRIPTION:
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
#  REQUIRED FILES:
#   Description:
#   Inpath:
#  OUTPUT:
#    Description: raw genotypes plus metadata for dartseq samples
#    Outpath:
#  DEPENDENCIES:
###.............................................................................
library(adegenet)
library(assertthat)
library(dplyr)
library(ape)

input_path <- paste0("data/intermediate/filt_genotypes.rds")
gen <- readRDS(file = input_path)

source("code/functions/bootstrap_trees.r")
source("code/functions/create_trees.r")

source("code/parameters/boot.r") #number of bootstrap

tr <- names(gen) %>% seq_along() %>%
sapply(function(j){
  if (adegenet::nLoc(gen[[j]]) < 100) type <- "usat"
    else if (adegenet::nLoc(gen[[j]]) > 100) type <- "snp"
  sp <- names(gen[j])
  temp <- list()
  original_trees <- create_trees(gen[[j]], type = type)
  boots_trees <- bootstrap_nj(gen[[j]], nboot, type = type)
  temp[[1]] <- c(original_trees, boots_trees)
  names(temp) <- sp
  temp
})


#calculate bootstrap support
bt <- names(tr) %>% seq_along() %>%
sapply(function(x){
  sp <- names(tr[x])
  temp <- list()
  phy <- reorder(tr[[x]][[1]], "postorder")
  ints <- phy$edge[, 2] > ape::Ntip(phy)
  ans <- ape::countBipartitions(phy, tr[[x]][-1])
  temp[[1]] <- c(nboot, ans[order(phy$edge[ints, 2])])
  names(temp) <- sp
  temp
})

#plot trees
  nodecol <- bt
for (i in seq_along(names(bt))){
  nodecol[[i]][nodecol[[i]] < tresh * nboot] <- "red"#color nodes bt BS
  nodecol[[i]][nodecol[[i]] != "red"] <- "black"
  pdf(file = paste0("data/intermediate/", names(bt)[[i]], "bs_tree.pdf"))
  plot(tr[[i]][[1]], "u", #plot fan
       use.edge.length = T,
       show.tip.label = T, cex = 0.4,
       edge.color = "black")
  nodelabels(pch = 20, col = nodecol[[i]], cex = 0.5) #issue assert that length of color nodes is the same as the number of nodes inf the tree
  add.scale.bar()
  title(paste0("NJ tree ", names(bt)[i], ", ", nLoc(gen[[i]]), " loci, ", nboot,
        " bootstrap")) #title plot
  dev.off()
      }
rm(i, tresh)

#save trees (first tree is the real tree)
saveRDS(tr, "data/intermediate/bt_trees.rds")
saveRDS(bt, "data/intermediate/boot_support.rds") #bootstrap support
