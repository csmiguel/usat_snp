###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: collect clumpak matrices
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(dplyr)
library(ggplot2)
library(cowplot)

source("code/parameters/plotting_par.r")

#1. create list with Major Cluster ancestries
qclumpak <-
  lapply(krange, function(k){
    distmatrix <-
    paste0("data/intermediate/clumpak_comparison/K=", k,
      "/CLUMPP.files/ClumppPairMatrix") %>%
    readLines() %>%
    textConnection %>%
    read.table
    runnames <- paste0("data/intermediate/clumpak_comparison/K=", k,
      "/CLUMPP.files/FilesToIndex") %>%
      readLines() %>%
      {sapply(strsplit(., "\t"), `[`, 2)}
    rownames(distmatrix) <- runnames
    colnames(distmatrix) <- runnames
    h <- reshape2::melt(as.matrix(distmatrix), varnames = c("row", "col"))
    hh <- h[as.numeric(h$row) > as.numeric(h$col), ]
    hh
    })
names(qclumpak) <- paste0("K", krange)
plot_data <-
  plyr::ldply(seq_along(qclumpak), function(i) {
    x <- qclumpak[[i]]
    x$cluster <- names(qclumpak)[i]
    x
  }) %>%
  as_tibble() %>%
    dplyr::mutate(comp = ifelse(grepl("usat", row) & grepl("usat", col), "usat",
                         ifelse(grepl("dart", row) & grepl("dart", col), "dart",
                          "both")))
#plots
labs_legend <- levels(as.factor(plot_data$comp)) %>%
  plyr::mapvalues(from = c("both", "usat", "dart"),
                  to = c("SNPs-microsatellites", "microsatellites", "SNPs"))

ggsave("data/final/SSC_markers.pdf", units = "cm", height = 12, width = 18)

ggplot(plot_data) +
    geom_point(aes(x = cluster, y = value, color = comp, shape = comp),
               alpha = 0.2, position = position_jitter(0.1), size = 4) +
    scale_color_manual(values = c("brown", colm),
                       guide = guide_legend(reverse = T),
                       name = "Pairwise comparison", labels = labs_legend) +
    scale_shape_manual(values = c(0, 1, 2), guide = guide_legend(reverse = T),
                       name = "Pairwise comparison", labels = labs_legend) +
    ylab("Symmetric Similarity Coefficient") + xlab(NULL) +
    theme(legend.position = c(0.3, 0.4),
    legend.box.background = element_rect(colour = "black", size = 0.2),
    legend.text = element_text(size = 10))
