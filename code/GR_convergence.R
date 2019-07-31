###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: cowplots convergence Gelman and Rubin's convergence diagnostic
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

library(dartR)
library(ggplot2)
library(flextable)
library(dplyr)
library(magrittr)

#load data
gen <- readRDS("data/intermediate/gen_consolidated_filtered.rds")
rhat_results <- readRDS("data/intermediate/rhat.rds")
source("code/functions/ft_helper.r")

#edit data
# remove upper_ci column and round numbers
rhat_results$dart %<>% dplyr::select(-upper_ci) %>%
  mutate(rhat = round(rhat, digits = 2)) %>%
  rename(dataset = id)
rhat_results$usat %<>% dplyr::select(-upper_ci) %>%
  mutate(rhat = round(rhat, digits = 2)) %>%
  rename(dataset = id)
# fix levels
levels(rhat_results$usat$dataset) <-
  plyr::mapvalues(
    x = levels(rhat_results$usat$dataset),
    from = c("usat_hyla", "usat_pelo"),
    to = c("H. molleri", "P. cultripes")
  )

##add number of loci to dart
h <- levels(rhat_results$dart$dataset)
h[grep("dart_hyla", h)] <- paste0("hyla_", adegenet::nLoc(gen$dart_hyla))
h[grep("dart_pelo", h)] <- paste0("pelo_", adegenet::nLoc(gen$dart_pelo))
levels(rhat_results$dart$dataset) <- h
rm(h, gen)

##Create tables
tableSNPs <-
  rhat_results$dart %>%
  dplyr::mutate(Species = stringr::str_extract(dataset, "[a-z]+")) %>%
  dplyr::mutate(Species = plyr::mapvalues(x = Species,
                                          from = c("hyla", "pelo"),
                                          to = c("H. molleri", "P. cultripes"))) %>%
  dplyr::mutate(dataset = stringr::str_extract(
    as.character(dataset), "[1-9].*") %>% as.numeric()) %>%
  {rhat_edit <<- dplyr::arrange(., dataset, k)} %>%
  reshape2::dcast(Species + dataset ~ k, value.var = "rhat") %>%
  dplyr::mutate( dataset = as.factor(dataset)) %>%
  dplyr::rename("No. of loci" = dataset) %>%
  flextable::flextable() %>%
  flextable::autofit() %>%
  flextable::colformat_num(digits = 2,
                           col_keys = as.character(unique(rhat_edit$k))) %>%
  flextable::align(j = 3:10, align = "center", part = "all") %>%
  flextable::merge_v(~Species) %>%
  flextable::hline(i = 7, border = officer::fp_border(color = "black"),
                   part = "body") %>%
  flextable::fix_border_issues() %>%
  flextable::italic(j = "Species") %>%
  {snpft <<- flextable::add_header_lines(., values = "A. SNPs")} %>%
  flextable::as_raster() %>%
  ft2ggplot()

#Table microsatellites
tableUSATs <-
  rhat_results$usat %>%
  {rhat_usat <<-dplyr::rename(., Species = dataset)} %>%
  reshape2::dcast(Species ~ k, value.var= "rhat") %>%
  flextable::flextable() %>%
  flextable::autofit() %>%
  flextable::colformat_num(digits = 2,
                           col_keys = as.character(unique(rhat_usat$k))) %>%
  flextable::italic(j = "Species") %>%
  flextable::align(j = 2:9, align = "center", part = "all") %>%
  {usatft <<- flextable::add_header_lines(., values = c("B. microsatellites"))} %>%
  flextable::as_raster() %>%
  ft2ggplot()
#create plot
p1 <-
  ggplot(data = rhat_edit, aes(x = k, y = rhat, color = as.factor(dataset))) +
  geom_line(size = 2) +
  ylab("Gelman and Rubin's convergence") +
  theme_classic() +
  xlab("K") +
  viridis::scale_color_viridis(name = "SNPs,\nNumber of loci",
                               discrete = T, option = "D") +
  facet_grid(~Species) +
  theme(strip.text.x = element_text(size = 10, face = "italic"))

## add line from microsatellites
p2 <-
  p1 +
  geom_line(data = rhat_usat, color = "grey", linetype = "solid", size = 1.5) +
  scale_linetype_manual(values = "1")

##
pdf(file = "data/final/GR_convergence.pdf", height = 9, width = 7)
cowplot::plot_grid(p2, tableSNPs, tableUSATs, nrow = 3,
                   rel_heights =  c(3, 3.5, 1), axis = "l")
dev.off()

#create export tables
usatft %<>% width(j = 2:9, 0.5)
snpft %<>% width(j = 3:10, 0.6)
officer::read_docx() %>%
  officer::body_add_par("Table Sx. Gelman and Rubin's convergence diagnostic") %>%
  body_add_flextable(value = snpft) %>%
  body_add_flextable(value = usatft) %>%
  print(target = "data/final/GR_convergence.docx")
