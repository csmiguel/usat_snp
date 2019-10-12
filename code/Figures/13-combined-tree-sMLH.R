###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: Plot trees with heterozygosity
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(ape)
library(dplyr)
library(adegenet)
library(grid)
library(ggplot2)
library(ggtree)
library(tidytree)

#bootstrap support
bt <- readRDS("data/intermediate/boot_support.rds")
#trees
tr <- readRDS("data/intermediate/bt_trees.rds")
#genotypes
gen <- readRDS("data/intermediate/gen_consolidated_filtered.rds")
#clumpak output: major cluster
qclumpak <- readRDS("data/intermediate/clumpak_major.rds")
#pop ids
meta <- readRDS("data/intermediate/table1.rds")
#heterozyosity
het <- readRDS("data/intermediate/sMLH.rds")

source("code/parameters/boot.r") #bootstrap parameters
source("code/functions/plot_nj_trees.r") #function for plotting
source("code/parameters/plotting_par.r") #plotting parameters

# 1. Create ggtrees with bootstrap support and structure ancestries.
hh <-
  #for each dataset do tree and facetted structure plots
  lapply(names(gen), function(z){
    tr1 <- tidytree::as.treedata(tr[[z]][[1]]) #convert tree to treedata
    het <- het[[z]]
    #compute ancestral states for heterozygosity
    fit <- phytools::fastAnc(tr1@phylo, het)
    #match heterozygosity to tips
    het2nodes <-
      dplyr::left_join(
        data.frame(tips = tr1@phylo$tip.label),
        data.frame(tips = names(het), het)) %>%
      dplyr::mutate(node = nodeid(tr1, tips)) %>%
      `rownames<-`(.$tips) %>%
      dplyr::select(node, het) %>%
      #rbind to fitted ancestries for internal nodes
      {rbind(., data.frame(node = names(fit), het = as.numeric(fit)))} %>%
      dplyr::mutate(node = as.numeric(node))
    assertthat::assert_that(all(!is.na(het2nodes)))
    tree <- full_join(tr1, het2nodes, by = "node")
    #create object color for nodes
    nodecol <- bt[[z]]
    coln[1] <- "white"
    nodecol[nodecol < thresh * nboot] <- coln[1]#color nodes bt BS
    nodecol[nodecol != coln[1]] <- coln[2]
    #create dataframe with pop ids (as in Table 1)
    ids <-
      data.frame(sample_id = tr[[z]][[1]]$tip.label) %>%
      dplyr::left_join(gen[[z]]$other$metadata %>%
                         dplyr::mutate(sample_id = as.character(sample_id )),
                       by = c("sample_id" = "sample_id")) %>%
      dplyr::left_join(meta, by = c("locality" = "locality")) %>%
      dplyr::select(sample_id, ID)
    #base tree
    p <-
      ggtree::ggtree(tree, color = "darkgrey") + #plot tree
      ggtree::geom_nodepoint(aes(color = het), size = 3) + #sMLH
      scale_color_gradientn(
        colours = rev(c("red", "orange", "green", "cyan", "blue")),
        breaks = seq(0, 2, 0.2)) +
      ggtree::geom_nodepoint(color = nodecol, size = 1) + #BS
      theme(legend.position = c(.05, .85)) +
      labs(color = "sMLH")
    # add tip labels
    p <- p %<+% ids + ggtree::geom_tiplab(aes(
      label = paste(ID, ids$sample_id, sep = "-")),
      size = 2, color = "black",
      align = T, linesize = 0.2, offset = 0)
    #facet plot with structure ancestries
    # for each K within each dataset
    for (k in 2:length(qclumpak[[z]])){
      #create matrix with ancestries
      hm <- qclumpak[[z]][[k]] %>%
        dplyr::mutate(tlab = adegenet::indNames(gen[[z]])) %>%
        reshape2::melt()
      #add structure facets to tree in a loop
      p <- ggtree::facet_plot(p, panel = paste0("K", k), data = hm,
                              geom = ggstance::geom_barh,
                              aes(x = value, fill = as.factor(variable)),
                              stat = "identity")
    }
    p2 <- p + theme_tree2(
      legend.position = c(.05, .9),
      legend.background = element_blank(),
      panel.spacing = unit(0.1, "lines"), #fix spacing between panels
            strip.background = element_rect(fill = NA, colour = NA), #remove bg
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()) + #x axis label size
      scale_x_continuous(expand = c(0, 0)) + #controls separation between panels
      scale_y_continuous(expand = c(0.01, 0.01)) + # gap bt top/bottom labs
      scale_fill_manual(values = str_col) + #color according to structure
      xlim_expand(c(max(p$data$x) * -0.1, max(p$data$x) * 1.3), "Tree") + #expand to fit indNames
      guides(fill = FALSE,
      colour = guide_colourbar(barheight = 5,
        draw.ulim = FALSE, draw.llim = FALSE))
    p2 <- ggtree::facet_labeller(p2,
                       c(Tree = plyr::mapvalues(#rename Tree with dataset
                         x = z,
                         from = names(gen),
                         to = c("SNPs H. molleri",
                                "SNPs P. cultripes",
                                "microsatellites H. molleri",
                                "microsatellites P. cultripes"))
                                 ))
    ggtree::facet_widths(p2, widths = c(2, rep(0.3, k - 1))) #fix panel widths
  })

# create multiplot
assertthat::assert_that(length(hh) == 4)
ggpubr::ggarrange(hh[[1]], hh[[2]], hh[[3]], hh[[4]],
                  ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave("data/final/combined_trees_heterozygosity.pdf", height = 16, width = 11)
