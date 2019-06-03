bt <- readRDS("data/intermediate/boot_support.rds")
tr <- readRDS("data/intermediate/bt_trees.rds")

#plot trees
  nodecol <- bt
for (i in seq_along(names(bt))){
  nodecol[[i]][nodecol[[i]] < tresh * nboot] <- "red"#color nodes bt BS
  nodecol[[i]][nodecol[[i]] != "red"] <- "black"
  pdf(file = paste0("data/final/", names(bt)[[i]], "bs_tree.pdf"))
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
