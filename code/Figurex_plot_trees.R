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
library(adegenet)
library(grid)
library(ggplot2)
library(ggtree)

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

source("code/parameters/boot.r") #bootstrap parameters
source("code/functions/plot_nj_trees.r") #function for plotting
source("code/parameters/plotting_par.r") #plotting parameters

# 1. Create ggtrees with bootstrap support and structure ancestries.

hh <-
  #for each dataset do tree and facetted structure plots
  lapply(names(gen), function(z){
  #create object color for nodes
  nodecol <- bt[[z]]
  nodecol[nodecol < thresh * nboot] <- coln[1]#color nodes bt BS
  nodecol[nodecol != coln[1]] <- coln[2]
  #create object with shape for nodes
  pchnod <- plyr::mapvalues(nodecol, coln, pchn) %>% as.numeric
  #create dataframe with pop ids (as in Table 1)
  ids <-
    data.frame(sample_id = tr[[z]][[1]]$tip.label) %>%
    dplyr::left_join(gen[[z]]$other$metadata %>%
                       dplyr::mutate(sample_id = as.character(sample_id )),
                     by = c("sample_id" = "sample_id")) %>%
    dplyr::left_join(meta, by = c("locality" = "locality")) %>%
    dplyr::select(sample_id, ID)
  #base tree
  p <- ggtree::ggtree(tr[[z]][[1]], color = "darkgrey", size = 0.2) + #plot tree
      ggtree::geom_nodepoint(color = nodecol, shape = pchnod, size = 1) #nodes
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
  p2 <- p + theme_tree2() +
    theme(panel.spacing = unit(0.1, "lines"), #fix spacing between panels
          strip.background = element_rect(fill = NA, colour = NA), #remove bg
          axis.text.x = element_text(size = 6)) + #x axis label size
    scale_x_continuous(#breaks = c(0, 0.25, 0.5, 0.75, 1), #manual ticks
                       #labels = c("","0.25", "", ".75", ""), #manual labels
                       expand = c(0, 0)) + #controls separation between panels
    scale_y_continuous(expand = c(0.01, 0.01)) + # gap bt top/bottom labs
    scale_fill_manual(values = str_col) + #color according to structure
    xlim_expand(c(0, max(p$data$x) * 1.3), "Tree") #expand tree to fit indNames
    p2 <- ggtree::facet_labeller(p2,
      c(Tree = plyr::mapvalues(#rename Tree with dataset
                  x = z,
                  from = names(gen),
                  to = c("SNPs H. molleri",
                         "SNPs P. cultripes",
                         "microtellites H. molleri",
                         "microtellites P. cultripes"))
                  ))
    ggtree::facet_widths(p2, widths = c(2, rep(0.3, k - 1))) #fix panel widths
  })

# create multiplot
assertthat::assert_that(length(hh) == 4)
ggpubr::ggarrange(hh[[1]], hh[[2]], hh[[3]], hh[[4]],
                  ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave("data/final/trees_structure.pdf", height = 16, width = 11)

# 2. Old trees with ape:
#rename datasets
names(tr) <-
  plyr::mapvalues(
    x = names(tr),
    from = names(tr),
    to = c("SNPs H. molleri", "SNPs P. cultripes",
    "microtellites H. molleri", "microtellites P. cultripes"))

#plot trees
plot_nj(bt = bt, thresh = thresh, coln = c("red", "black"),
  tips = "none", typeP = "unrooted", K = "K2", output = ".pdf")
#plots using ancestries from major cluster in K = 2
plot_nj(bt = bt, thresh = thresh, coln = c("red", "black"),
  tips = "ancestry", typeP = "phylogram", K = "K2", output = ".pdf")
#plots using ancestries from major cluster in K = 4
plot_nj(bt = bt, thresh = thresh, coln = c("red", "black"),
  tips = "ancestry", typeP = "phylogram", K = "K4", output = ".pdf")
