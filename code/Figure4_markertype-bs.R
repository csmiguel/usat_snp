###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: Plot results from bootstrap analysis
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

library(dplyr)
library(ggplot2)
library(cowplot)

bstree <- readRDS("data/intermediate/bs_treemetrics.rds")
source("code/parameters/plotting_par.r")
source("code/functions/plot_bs.r")

#rename datasets
bstree
names(bstree) <-
  plyr::mapvalues(
    x = names(bstree),
    from = names(bstree),
    to = c("SNPs H. molleri", "SNPs P. cultripes",
           "microtellites H. molleri", "microtellites P. cultripes"))
#1. Old plots BW

# plot parameters
x_lab <- c("Distance to parent node", "Distance from root")
max_dnr <- sapply(bstree, function(x) x$dnr) %>% unlist %>% max
max_dnp <- sapply(bstree, function(x) x$dnp) %>% unlist %>% max
x_lim <- matrix(c(0, max_dnp, 0, max_dnr), ncol = 2)
#loop plots
for (i in 1:2){
  pdf(file = paste0("data/final/",
                    gsub(" ", "_", x_lab[i]), ".pdf"),
      height = 7, width = 8.8)
  par(mfrow = c(2, 2), las = 1)
  for (p in 1:length(bstree)){
    hh <- bstree[[p]]
    plot(as.numeric(as.matrix(hh[, i])),
         hh$bs_tm, pch = 20,
         xlab = x_lab[i],
         ylab = "Bootstrap support",
         main = names(bstree)[p],
         xlim = x_lim[, i],
         ylim = c(0, 1))
  }
  dev.off()
}

#2. ggplot2 plots (facetted)
# format data for ggplot2
plot_data <-
  seq_along(bstree) %>%
  lapply(function(x){
    bstree[[x]] %>% mutate(dataset = as.factor(names(bstree)[x]))
  }) %>% {do.call(rbind, .)} %>%
  dplyr::mutate(species = as.character(dataset) %>%
                  stringr::str_remove("[aA-zA]++ ")) %>%
  dplyr::mutate(marker = as.character(dataset) %>%
                  stringr::str_remove(" .*$")) %>%
  {reshape2::melt(., id.vars = names(.)[-grep("dnr|dnp", names(.))])} %>%
  tibble::as_tibble() %>%
  dplyr::select(-dataset) %>%
  dplyr::rename(dist = variable)

#labels distances
labels_distance <- c(dnp = "Dnp", dnr = "Dnr")

#produce plots
plot_bs(mode = "bw")
plot_bs(mode = "color")
