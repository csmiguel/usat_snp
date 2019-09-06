###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: Plot results from bootstrap analysis per marker
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
