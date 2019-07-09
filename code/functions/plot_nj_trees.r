###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: plot nj trees
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

plot_nj <-  function(bt = bt, thresh = thresh, coln = c("red", "black"),
                tips = c("names", "ancestry"), typeP = c("phylogram", "unrooted"),
                output = "bs_tree.pdf", pchn = c(20, 18), K = "K2"){
  #bt is vector with bootstrap supports
  #thresh is the threshold to color nodes
  #coln are the colors for below, above threshold
  #tips, names: show individual labels, ancestry: shows pie with ancestry.
  nodecol <- bt
  if(tips == "names") tip <- T else tip <- F
  for (i in seq_along(bt)){
    nodecol[[i]][nodecol[[i]] < thresh * nboot] <- coln[1]#color nodes bt BS
    nodecol[[i]][nodecol[[i]] != coln[1]] <- coln[2]
    pchnod <- plyr::mapvalues(nodecol[[i]], coln, pchn) %>% as.numeric
    pdf(file = paste0("data/final/",
      paste(names(bt)[[i]], K, typeP, output, sep = "_")))
    plot(tr[[i]][[1]], type = typeP,
         use.edge.length = T,
         show.tip.label = tip, cex = 0.4,
         edge.color = "grey")
    #issue assert that length of color nodes is the same as the number of nodes inf the tree
    ape::nodelabels(pch = pchnod, col = nodecol[[i]], cex = 0.8)
    add.scale.bar()
    title(paste0("NJ tree ", names(tr)[i], ", ", adegenet::nLoc(gen[[i]]),
                 " loci, ", nboot, " bootstrap"), cex.main = 0.7) #title plot
    if (tips == "ancestry"){
      #get ancestries from clumpak object
      h <- qclumpak[[names(gen)[i]]][[K]] %>%
        as.matrix %>%
        `rownames<-`(adegenet::indNames(gen[[names(gen)[i]]]))
      assertthat::assert_that(all(tr[[i]][[1]]$tip.label %in% rownames(h)))
      #sort ancestries according to tree labels
      h <- h[match(tr[[i]][[1]]$tip.label, rownames(h)), ]
      #vector with pop ids
      ids <-
        h %>% as.data.frame() %>% tibble::rownames_to_column() %>%
        dplyr::select(1) %>%
        dplyr::left_join(gen[[names(gen)[i]]]$other$metadata %>%
        dplyr::mutate(sample_id = as.character(sample_id )),
          by = c("rowname" = "sample_id")) %>%
        dplyr::left_join(meta, by = c("locality" = "locality")) %>%
        dplyr::select(ID) %>% unlist
      #plot tip labels
        #create offset which is a proportion from the plot window
      offset1 <- get("last_plot.phylo", envir = .PlotPhyloEnv)$x.lim[2] * 0.02
      tiplabels(pie = h, cex = 0.3, piecol = topo.colors(ncol(h)), offset = 0)
      tiplabels(cex = 0.4, text = ids, frame = "none", offset = offset1)
    }
    dev.off()
      }
    }
