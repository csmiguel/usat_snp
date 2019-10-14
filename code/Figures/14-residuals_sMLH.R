###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: Plot trees residuals sMLH
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
het <- readRDS("data/intermediate/sMLH.rds") %>%
  {.[grep("pelo", names(.))]}

#residuals sMLH
res <- readRDS("data/intermediate/residuals-sMLH.rds")


source("code/parameters/boot.r") #bootstrap parameters
source("code/functions/plot_nj_trees.r") #function for plotting
source("code/parameters/plotting_par.r") #plotting parameters

# 1. Create ggtrees with bootstrap support and structure ancestries.

tr1 <- tidytree::as.treedata(tr[["dart_pelo"]][[1]]) #convert tree to treedata
#compute ancestral states for residual heterozygosity
names(res) <- tr1@phylo$tip.label
fit <- phytools::fastAnc(tr1@phylo, res)
#match residuals herozygosity to tips
res2nodes <-
  dplyr::left_join(
    data.frame(tips = tr1@phylo$tip.label),
    data.frame(tips = names(res), res)) %>%
  dplyr::mutate(node = nodeid(tr1, tips)) %>%
  `rownames<-`(.$tips) %>%
  dplyr::select(node, res) %>%
  #rbind to fitted ancestries for internal nodes
  {rbind(., data.frame(node = names(fit), res = as.numeric(fit)))} %>%
  dplyr::mutate(node = as.numeric(node))
assertthat::assert_that(all(!is.na(res2nodes)))
tree <- full_join(tr1, res2nodes, by = "node")
#create dataframe with pop ids (as in Table 1)
ids <-
  data.frame(sample_id = tr[["dart_pelo"]][[1]]$tip.label) %>%
  dplyr::left_join(gen[["dart_pelo"]]$other$metadata %>%
                     dplyr::mutate(sample_id = as.character(sample_id )),
                   by = c("sample_id" = "sample_id")) %>%
  dplyr::left_join(meta, by = c("locality" = "locality")) %>%
  dplyr::select(sample_id, ID)
#base tree
p <-
  ggtree::ggtree(tree, color = "darkgrey") + #plot tree
  ggtree::geom_nodepoint(aes(color = res), size = 3) + #sMLH
  ggtree::geom_nodepoint(size = 3, shape = 1) +
  scale_color_gradientn(
    colours = c("blue", "white", "red")) +
  labs(color = "residuals sMLH \nSNPs-microsatellites")
# add tip labels
p <- p %<+% ids + ggtree::geom_tiplab(aes(
  label = paste(ID, ids$sample_id, sep = "-")),
  size = 2, color = "black",
  align = T, linesize = 0.2, offset = 0)
#facet plot with structure ancestries
# for each K within each dataset
for (k in 2:length(qclumpak[["dart_pelo"]])){
  #create matrix with ancestries
  hm <- qclumpak[["dart_pelo"]][[k]] %>%
    dplyr::mutate(tlab = adegenet::indNames(gen[["dart_pelo"]])) %>%
    reshape2::melt()
  #add structure facets to tree in a loop
  p <- ggtree::facet_plot(p, panel = paste0("K", k), data = hm,
                          geom = ggstance::geom_barh,
                          aes(x = value, fill = as.factor(variable)),
                          stat = "identity")
}
p2 <- p + theme_tree2(
  legend.position = c(.08, .92),
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
         colour = guide_colourbar(barheight = 5))
p2 <- ggtree::facet_labeller(p2,
                             c(Tree = plyr::mapvalues(#rename Tree with dataset
                               x = "dart_pelo",
                               from = names(gen),
                               to = c("SNPs H. molleri",
                                      "SNPs P. cultripes",
                                      "microsatellites H. molleri",
                                      "microsatellites P. cultripes"))
                             ))
p3 <- ggtree::facet_widths(p2, widths = c(2, rep(0.3, k - 1))) #fix panel widths

#plot Correlation
assertthat::assert_that(
  all(names(het[["usat_pelo"]]) == names(het[["dart_pelo"]])))

plot_het <-
  seq_along(het) %>%
  sapply(function(x){
    het[[x]]
  }) %>% as.data.frame
names(plot_het) <- names(het)

p4 <-
  ggplot(plot_het, aes(usat_pelo, dart_pelo)) +
  geom_point() +
  labs(x = "sMLH microsatellites", y = "sMLH SNPs") +
  stat_smooth(method="lm", formula = y ~ x, colour = "lightblue") +
  theme_classic()

# multiplot
top_row <- cowplot::plot_grid(p4, NULL, rel_widths = c(1, 1))
cowplot::plot_grid(top_row, p3,
                   labels = "AUTO",
                   ncol = 1,
                   rel_heights = c(1, 3))
ggsave("data/final/residuals_sMLH.pdf", height = 16, width = 11)
