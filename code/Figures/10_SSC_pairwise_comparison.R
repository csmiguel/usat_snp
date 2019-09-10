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

#import data with Symmetric similarity coefficients (SSC)
plot_data <- readRDS("data/intermediate/SSC.rds")

source("code/parameters/plotting_par.r")

#create labels from legend
labs_legend <- levels(as.factor(plot_data$comp)) %>%
  plyr::mapvalues(from = c("both", "usat", "dart"),
                  to = c("SNPs-microsatellites",
                  "microsatellites-microsatellites",
                  "SNPs-SNPs"))
#color plot
ggplot(plot_data) +
    geom_point(aes(x = cluster, y = value, color = comp, shape = comp),
               alpha = 1, position = position_jitter(0.1), size = 4) +
    scale_color_manual(values = c("brown", colm),
                       guide = guide_legend(reverse = T),
                       name = "Pairwise comparison", labels = labs_legend) +
    scale_shape_manual(values = c(0:2), guide = guide_legend(reverse = T),
                       name = "Pairwise comparison", labels = labs_legend) +
    ylab("Symmetric Similarity Coefficient") + xlab(NULL) +
    theme(legend.position = c(0.3, 0.4),
    legend.box.background = element_rect(colour = "black", size = 0.2),
    legend.text = element_text(size = 10))

ggsave("data/final/SSC_markers_color.pdf", units = "cm", height = 12, width = 18)

#BW plot
# create colors for point filling
colbw1 <- c(grey(0.7),  "white", grey(0.3))

ggplot(plot_data) +
    geom_point(aes(x = cluster, y = value, color = comp, shape = comp, fill = comp),
               alpha = 1, position = position_jitter(0.1), size = 4) +
    scale_color_manual(values = rep("black", 3),
                       guide = guide_legend(reverse = T),
                       name = "Pairwise comparison", labels = labs_legend) +
    scale_fill_manual(values = colbw1,
                      guide = guide_legend(reverse = T),
                      name = "Pairwise comparison", labels = labs_legend) +
    scale_shape_manual(values = c(22, 21, 24), guide = guide_legend(reverse = T),
                       name = "Pairwise comparison", labels = labs_legend) +
    ylab("Symmetric Similarity Coefficient") + xlab(NULL) +
    theme(legend.position = c(0.3, 0.4),
          legend.box.background = element_rect(colour = "black", size = 0.2),
          legend.text = element_text(size = 10))

ggsave("data/final/SSC_markers_BW.pdf", units = "cm", height = 12, width = 18)
