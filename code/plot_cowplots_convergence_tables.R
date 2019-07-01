###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: cowplots convergence Gelman and Rubin's convergence diagnostic
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

library(magrittr)
library(flextable)
library(grid)
library(dplyr)
library(ggplot2)


rhat_results <- readRDS("data/intermediate/rhat.rds")

#remove upper_ci column and round numbers
rhat_results$dart %<>% dplyr::select(-upper_ci) %>%
  mutate(rhat = round(rhat, digits = 2)) %>%
  rename(dataset = id)
rhat_results$usat %<>% dplyr::select(-upper_ci) %>%
  mutate(rhat = round(rhat, digits = 2)) %>%
  rename(dataset = id)


#replacements in dataset names
levels(rhat_results$dart$dataset) <-
  plyr::mapvalues(
  x = levels(rhat_results$dart$dataset),
  from = c("1000", "10000", "200", "20000", "3000", "500", "5000",
  "dart_hyla", "dart_pelo"),
  to = c("subset 1000", "subset 10000", "subset 200", "subset 20000",
  "subset 3000", "subset 500", "subset 5000", "H. molleri", "P. cultripes")
  )

levels(rhat_results$usat$dataset) <-
  plyr::mapvalues(
  x = levels(rhat_results$usat$dataset),
  from = c("usat_hyla", "usat_pelo"),
  to = c("H. molleri", "P. cultripes")
  )

#reorder levels
rhat_results$dart$dataset <-
factor(rhat_results$dart$dataset,
   levels = c("H. molleri", "P. cultripes", "subset 200", "subset 500",
   "subset 1000", "subset 3000", "subset 5000", "subset 10000", "subset 20000"))

#reshape tables
rusat <-
  rhat_results$usat %>%
  reshape2::dcast(dataset ~ k)
rdart <-
  rhat_results$dart %>%
  reshape2::dcast(dataset ~ k)

#create flextables
ft_raster_dart <-
  flextable(rdart) %>%
  autofit() %>%
  add_header_row(values = c("", "K"), colwidths = c(5, 4)) %>%
  as_raster()

ft_raster_usat <-
  flextable(rusat) %>%
  autofit() %>%
  add_header_row(values = c("", "K"), colwidths = c(5, 4)) %>%
  as_raster()

#ggplots for flextables
pdart <- ggplot() +
  theme_void() +
  annotation_custom(grid::rasterGrob(ft_raster_dart),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
pusat <- ggplot() +
  theme_void() +
  annotation_custom(grid::rasterGrob(ft_raster_usat),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

#plot lines Gelman and Rubin's convergence
p_usat2 <-
  ggplot(data = rhat_results$usat, aes(x = k, y = rhat, color = dataset)) +
  geom_line() +
  ggtitle("microsatellites") +
  ylab("Gelman and Rubin's convergence")
p_dart2 <-
  ggplot(data = rhat_results[[2]], aes(x = k, y = rhat, color = dataset)) +
  geom_line() +
  ggtitle("SNPs") +
  ylab("Gelman and Rubin's convergence")
##
#save plots to files
# usats
pdf(paste0("data/final/GR_convergence_usats.pdf"), height = 6, width = 8)
cowplot::plot_grid(p_usat2, pusat, nrow = 2, ncol = 1, rel_heights = c(1, 1))
dev.off()

#snps
pdf(paste0("data/final/GR_convergence_snps.pdf"), height = 6, width = 8)
cowplot::plot_grid(p_dart2, pdart, nrow = 2, ncol = 1, rel_heights = c(1, 1))
dev.off()
