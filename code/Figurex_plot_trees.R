###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: Plot results from bootstrap analysis
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(ape)
library(dplyr)

#bootstrap support
bt <- readRDS("data/intermediate/boot_support.rds")
#trees
tr <- readRDS("data/intermediate/bt_trees.rds")
#genotypes
gen <- readRDS("data/intermediate/gen_consolidated_filtered.rds")
#best files as determined by KFinder
best_files <- readRDS("data/intermediate/best_files_kfinder.rds")

source("code/parameters/boot.r")

#rename datasets
names(tr) <-
  plyr::mapvalues(
    x = names(tr),
    from = names(tr),
    to = c("SNPs H. molleri", "SNPs P. cultripes",
    "microtellites H. molleri", "microtellites P. cultripes"))

#plot trees
plot_nj <-  function(bt = bt, thresh = thresh, coln = c("red", "black"),
                tips = c("names", "ancestry"), output = "bs_tree.pdf"){
  #bt is vector with bootstrap supports
  #thresh is the threshold to color nodes
  #coln are the colors for below, above threshold
  #tips, names: show individual labels, ancestry: shows pie with ancestry.
  nodecol <- bt
  if(tips == "names") tip <- T else tip <- F
  for (i in seq_along(bt)){
    nodecol[[i]][nodecol[[i]] < thresh * nboot] <- coln[1]#color nodes bt BS
    nodecol[[i]][nodecol[[i]] != coln[1]] <- coln[2]
    pdf(file = paste0("data/final/", names(bt)[[i]], output))
    plot(tr[[i]][[1]], "u", #plot fan
         use.edge.length = T,
         show.tip.label = tip, cex = 0.4,
         edge.color = "grey")
    #issue assert that length of color nodes is the same as the number of nodes inf the tree
    nodelabels(pch = 20, col = nodecol[[i]], cex = 0.5)
    if (tips == "ancestry"){
      h <- pophelper::readQStructure(best_files[8], indlabfromfile = T)[[1]] %>%
       as.matrix
      assertthat::assert_that(all(tr[[1]][[1]]$tip.label %in% rownames(h)))
      h <- h[match(tr[[1]][[1]]$tip.label, rownames(h)), ]
      tiplabels(pie = h, cex = 0.3, piecol = topo.colors(ncol(h)))
    }
    add.scale.bar()
    title(paste0("NJ tree ", names(tr)[i], ", ", adegenet::nLoc(gen[[i]]),
    " loci, ", nboot, " bootstrap"), cex.main = 0.7) #title plot
    dev.off()
      }
    }
