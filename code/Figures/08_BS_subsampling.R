# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: Plot results from subsampling SNPs from Hyla and Pelobates
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

library(dplyr)
library(magrittr)
library(dartR)
library(ggplot2)

gen <- readRDS("data/intermediate/gen_consolidated_filtered.rds")
source("code/parameters/plotting_par.r")

nloc_pelo <- gen$dart_pelo %>% adegenet::nLoc()
nloc_hyla <- gen$dart_hyla %>% adegenet::nLoc()
nloc_mhyla <- gen$usat_hyla %>% adegenet::nLoc()
nloc_mpelo <- gen$usat_pelo %>% adegenet::nLoc()

#create friendly ggplot2 dataframe with bootstrap data for all datasets
bt <-
  #read subsampled datasets
  readRDS("data/intermediate/s_bs_treemetrics.rds") %>%
  #concatenate lists within each dataset(all, subhyla, subpelo)
  {c(list("all" = readRDS("data/intermediate/bs_treemetrics.rds")), .)} %>%
  reshape2::melt() %>%
  dplyr::filter(variable == "bs_tm") %>%
  dplyr::as_tibble() %>%
  dplyr::select(-1) %>%
  dplyr::mutate(species = c("pelo", "hyla")[(grepl("hyla",
                                            L2) | grepl("hyla", L1)) + 1]) %>%
  dplyr::mutate(marker = c("snp", "usat")[(grepl("usat", L2)) + 1]) %>%
  dplyr::mutate(
    nloci = stringr::str_replace(L2, "dart_hyla", as.character(nloc_hyla)) %>%
      stringr::str_replace("dart_pelo", as.character(nloc_pelo)) %>%
      stringr::str_replace("usat_pelo", as.character(nloc_mpelo)) %>%
      stringr::str_replace("usat_hyla", as.character(nloc_mhyla)) %>%
      as.numeric()
  ) %>%
  dplyr::select(-L1, -L2)

xlabss <- as.numeric(c("50", "200", "500", "1000", "3000", "5000", "10000"))
assertthat::assert_that(all(xlabss %in% bt$nloci))

##ggplot2
plot_bs_loci_subsampling <- function(mode = c("bw", "color")){
  #mode is color or black and white
  if (mode == "color"){
    colbx <- c("#347b34", "#c23803")
    colpt <- colsp
  } else if (mode == "bw"){
    colbx <- c("black", "darkgrey")
    colpt <- colmbw
  }

p1 <-
  bt %>% dplyr::filter(marker == "snp") %>%
    dplyr::mutate(nloci = as.factor(nloci)) %>%
  ggplot(aes(x = as.numeric(as.character(nloci)), y = value, fill = species)) +
    geom_point(data = . %>% filter(species == "hyla"),
     aes(x = as.numeric(as.character(nloci)) * .9,
     y = value, colour = species), position = position_jitter(width = .03),
     size = 1, shape = 16, alpha = .5) +
    geom_point(data = . %>% filter(species == "pelo"),
     aes(x = as.numeric(as.character(nloci)) * 1.1,
     y = value, colour = species), position = position_jitter(width = .03),
     size = 1, shape = 16, alpha = .5) +
  geom_boxplot(data = . %>% filter(species == "hyla"),
     aes(group = nloci, x = as.numeric(as.character(nloci)) * .9, y = value),
     outlier.shape = NA, alpha = .5, colour = colbx[1], width = .05, coef = 0) +
  geom_boxplot(data = . %>% filter(species == "pelo"),
     aes(group = nloci, x = as.numeric(as.character(nloci)) * 1.1, y = value),
     colour = colbx[2], outlier.shape = NA, alpha = .5, width = .05, coef = 0) +
  scale_color_manual(values = colpt, name = "Species",
     labels = c("H. molleri", "P. cultripes")) +
  scale_fill_manual(values = colpt) +
  stat_smooth(method = "lm", formula = y ~ log(x),
    colour = "black", size = 0.5) +
  theme_classic(base_size = 8) +
  ylab("Bootstrap support") +
    scale_x_log10("natural logarithm of no. of SNPs",
              breaks = xlabss, labels = xlabss) +
  scale_y_continuous(expand = c(0.02, 0)) +
  guides(fill = FALSE,
         color = guide_legend(
           label.theme = element_text(face = "italic", size = 7))) +
  theme(legend.position = c(0.83, 0.12),
        legend.spacing.y = unit(0.1, "line"),
        legend.key.height = unit(0.5, "line"),
        legend.background = element_rect(fill = "transparent"))

#Plot microsatellite data
p2 <-
bt %>% dplyr::filter(marker == "usat") %>%
  dplyr::mutate(nloci = as.factor(nloci)) %>%
  ggplot(aes(x = as.numeric(as.character(nloci)), y = value, fill = species)) +
  geom_point(aes(x = species, y = value, colour = species),
             position = position_jitter(width = .1),
             size = 1, shape = 16, alpha = .5) +
  geom_boxplot(aes(x = species, y = value, colour = species),
             outlier.shape = NA, alpha = .5, width = .2,
             colour = colbx, coef = 0) +
  theme_classic(base_size = 8) +
  scale_color_manual(values = colpt) +
  scale_fill_manual(values = colpt) +
  scale_x_discrete(name = "Microsatellites",
    labels = c("H. molleri", "P. cultripes")) +
  scale_y_continuous(expand = c(0.02, 0)) +
  guides(fill = FALSE,
         color = F) +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(face = "italic"))
cowplot::plot_grid(p1, NULL, p2, nrow = 1,
  align = "h", rel_widths = c(1, -.05, .3))
}

#generate plots
plot_bs_loci_subsampling(mode = "bw")
ggsave("data/final/BS_loci_subsampling_bw.pdf",
  width = 120, height = 70, units = "mm")

plot_bs_loci_subsampling(mode = "color")
ggsave("data/final/BS_loci_subsampling_color.pdf",
       width = 120, height = 70, units = "mm")
