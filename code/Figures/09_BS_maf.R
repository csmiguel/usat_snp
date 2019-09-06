##
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: Plot results from MAF from Hyla and Pelobates
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

library(dplyr)
library(ggplot2)

bstree <- readRDS("data/intermediate/maf_bs_treemetrics.rds")

source("code/parameters/plotting_par.r")

#create friendly ggplot2 dataframe with bootstrap data for all datasets
bt <-
  #read subsampled datasets
  readRDS("data/intermediate/maf_bs_treemetrics.rds") %>%
  reshape2::melt() %>%
  dplyr::filter(variable == "bs_tm") %>%
  dplyr::as_tibble() %>%
  dplyr::select(-1) %>%
  dplyr::rename(maf = L2) %>%
  dplyr::mutate(species = c("pelo", "hyla")[(grepl("hyla", L1)) + 1]) %>%
  dplyr::select(-L1)

##ggplot2
plot_bs_maf <- function(mode = c("bw", "color")){
  jitbox <- 0.003
  #mode is color or black and white
  if (mode == "color"){
    colbx <- c("#347b34", "#c23803")
    colpt <- colsp
  } else if (mode == "bw"){
    colbx <- c("black", "darkgrey")
    colpt <- colmbw
  }

  xlabss <- unique(bt$maf) %>% as.numeric() %>% sort()

    bt %>% dplyr::mutate(maf = as.factor(maf)) %>%
    ggplot(aes(x = as.numeric(as.character(maf)), y = value, fill = species)) +
    geom_point(data = . %>% filter(species == "hyla"),
       aes(x = as.numeric(as.character(maf)) - jitbox,
       y = value, colour = species), position = position_jitter(width = .001),
       size = 1, shape = 16, alpha = .5) +
    geom_point(data= . %>% filter(species == "pelo"),
       aes(x = as.numeric(as.character(maf)) + jitbox,
       y = value, colour = species), position = position_jitter(width = .001),
       size = 1, shape = 16, alpha = .5) +
    geom_boxplot(data = . %>% filter(species == "hyla"),
       aes(group = maf, x = as.numeric(as.character(maf)) - jitbox, y = value),
       outlier.shape = NA, alpha = .5, colour = colbx[1],
       width = .004, coef = 0) +
    geom_boxplot(data = . %>% filter(species == "pelo"),
       aes(group = maf, x = as.numeric(as.character(maf)) + jitbox, y = value),
       colour = colbx[2], outlier.shape = NA, alpha = .5,
       width = .004, coef = 0) +
    scale_color_manual(values = colpt, name = "Species",
                       labels = c("H. molleri", "P. cultripes")) +
    scale_fill_manual(values = colpt) +
    stat_smooth(method = "lm", formula = y ~ x, ##
                colour = "black", size = 0.5) +
    theme_classic(base_size = 8) +
    ylab("Bootstrap support") +
    scale_x_continuous("Minor Allele Frequency",
                  breaks = xlabss, labels = xlabss) +
    scale_y_continuous(expand = c(0.02, 0)) +
    guides(fill = FALSE,
           color = guide_legend(
             label.theme = element_text(face = "italic", size = 7))) +
    theme(legend.position = "right",
           legend.spacing.y = unit(0.2, "line"),
           legend.key.height = unit(0.7, "line"),
           legend.background = element_rect(fill = "transparent"))
}

#generate plots
plot_bs_maf(mode = "bw")
ggsave("data/final/BS_maf_bw.pdf", width = 100, height = 56, units = "mm")

plot_bs_maf(mode = "color")
ggsave("data/final/BS_maf_color.pdf", width = 100, height = 56, units = "mm")
