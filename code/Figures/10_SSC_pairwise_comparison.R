###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: Plot distances from Pairwise comparisons
# see Jakobsson and Rosenberg 2002 10.1093/bioinformatics/btm233
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

library(ggplot2)
library(cowplot)
library(dplyr)

#import data with Symmetric similarity coefficients (SSC)
ssc_pelo <- readRDS("data/intermediate/SSC_pelo.rds") %>%
  {.$species <- "pelo";.}
ssc_hyla <- #reshape data to have the same format as ssc_pelo
  readRDS("data/intermediate/ssc_hyla.rds") %>%
  reshape2::melt() %>%
  dplyr::select(1, 2, 5, 6, 3) %>%
  dplyr::as_tibble() %>%
  setNames(names(ssc_pelo)[1:ncol(ssc_pelo) -1]) %>%
  {.$species <- "hyla";.}
plot_data <- rbind(ssc_pelo, ssc_hyla) %>%
              dplyr::mutate(species = as.factor(species))

source("code/parameters/plotting_par.r")

#create labels from legend

labs_legend <- levels(as.factor(plot_data$comp)) %>%
  plyr::mapvalues(from = c("both", "usat", "dart"),
                  to = c("SNPs-microsatellites",
                  "microsatellites-microsatellites",
                  "SNPs-SNPs"))
sp_labs <- c("Hyla molleri", "Pelobates cultripes") %>%
  setNames(levels(plot_data$species))
cowplot::theme_cowplot()

#plot
p_vertical <-
  ggplot(plot_data) +
  geom_point(aes(x = reorder(cluster, desc(cluster)), y = value, color = comp,
                 shape = comp),
             alpha = 1, position = position_jitter(0.1), size = 2) +
  scale_color_manual(values = c("brown", colm),
                     guide = guide_legend(reverse = T),
                     name = "Pairwise comparison", labels = labs_legend) +
  scale_shape_manual(values = c(0:2), guide = guide_legend(reverse = T),
                     name = "Pairwise comparison", labels = labs_legend) +
  ylab("Symmetric Similarity Coefficient") + xlab(NULL) +
  coord_flip() +
  facet_wrap(~species, dir = "v",
             labeller = labeller(species = sp_labs)) +
  theme_classic() +
  theme(legend.position = c(0.3, 0.4),
        legend.title = element_blank(),
        axis.text = element_text(colour = "black"),
        text = element_text(size = 10),
        strip.text.x = element_text(face = "italic"),
        strip.text.y = element_text(face = "italic"),
        strip.background = element_rect(
          color = NA, fill = NA, size = 1.5, linetype = "solid"),
        legend.key.size = unit(0, "lines"))
p_vertical
leg <- cowplot::get_legend(p_vertical)

p_vertical <- p_vertical + theme(legend.position = "none")
cowplot::plot_grid(ggdraw(leg), p_vertical, ncol = 1, rel_heights = c(0.1, 1))
ggsave("data/final/SSC_markers_color_flipped.pdf", units = "in",
       height = 7.5, width = 2.5)
