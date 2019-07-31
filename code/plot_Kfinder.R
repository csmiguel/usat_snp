###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: plot
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(dplyr)
library(magrittr)

#load kfinder results
kfinder <- readRDS("data/intermediate/kfinder_tidy.rds")
evanno <- readRDS("data/intermediate/evanno_tidy.rds")
source("code/functions/plot_kfinder.r")
source("code/parameters/plotting_par.r")
source("code/functions/ft_helper.r")

# 1. create word table with data
evanno %<>%
  dplyr::filter(variable == "L(K)") %>%
  dplyr::filter(!grepl("[0-9]", dataset)) %>%
  dplyr::mutate(toremove = paste0(dataset, K))

kfinder %>%
  dplyr::filter(!grepl("[0-9]", dataset)) %>%
  dplyr::mutate(toremove = paste0(dataset, K)) %>%
  dplyr::left_join(evanno, by = c("toremove" = "toremove")) %>%
  dplyr::select(species.x, marker.x, K.x, Pritchard, sd, Evanno, Parsimony) %>%
  {bestKs <<- data.table::setnames(.,
   c("species.x", "marker.x", "K.x", "sd", "Pritchard", "Evanno", "Parsimony"),
   c("species", "marker", "K", "sd Pr[X|K]",
     "Pr[X|K]", "Evanno ΔK", "Parsimony index")); bestKs} %>%
  reshape2::melt(id.vars = c("species", "marker", "K")) %>%
  reshape2::dcast(species + marker + variable ~ K) %>%
  dplyr::mutate_each_(funs(as.numeric), 4:11) %>%
  {bestk_reshaped <<- dplyr::rename(., metric = variable)} %>%
  flextable::flextable() %>%
  flextable::merge_v(~ species + marker) %>%
  flextable::fix_border_issues() %>%
  flextable::hline(i = 8, part = "body",
                   border = officer::fp_border(color = "black", width = 1)) %>%
  flextable::hline(i = c(4, 12), part = "body",
                   border = officer::fp_border(color = "black", width = 1)) %>%
  flextable::italic(j = 1) %>%
  flextable::colformat_num(col_keys = as.character(1:8),
                           digits = 2, big.mark = "") %>%
  flextable::colformat_num(i = c(1, 5, 9, 13), col_keys = as.character(1:8),
                           digits = 0, big.mark = "") %>%
  flextable::autofit() %>%
  flextable::width(j = 4:11, 0.9) %>%
  flextable::width(j = 1, 0.9) %>%
  flextable::width(j = 2, 1.2) %>%
  ft2word(path = "data/raw/Docxlandscape.docx",
          name = "data/final/bestK_metrics.docx")

#2. Plot bestKs

library(ggplot2)
pdf.options(encoding='ISOLatin2.enc')
ggplot_bestk <- function(colo = NULL, sizep = 7, aspp = 1.6){
  #colo must be colmbw or colm
  toplot <- reshape2::melt(bestk_reshaped) %>%
    dplyr::mutate(K = as.numeric(as.character(variable))) %>%
    dplyr::select(-variable) %>%
    dplyr::filter(metric == "Evanno ΔK" | metric == "Parsimony index")

  ggplot(data = toplot, aes(x = K, y = value)) +
    geom_line(aes(linetype = marker, color = marker)) +
    geom_point(aes(shape = marker, color = marker)) +
    theme_classic() +
    theme(strip.text.x = element_text(face = "italic"),
          strip.background = element_rect(color = "white"),
          panel.grid.major.x = element_line(
            size = .1, color = "grey", linetype = "dotted")) +
    scale_linetype_manual(values = ltym) +
    scale_color_manual(values = colo) +
    scale_shape_manual(values = pchm) +
    facet_grid(metric ~ species, scales = "free_y")
  ggsave(paste0("data/final/Evanno_sd_", substitute(colo),
                ".pdf"), width = sizep, height = sizep / aspp)

  ggplot(data = bestKs, aes(x = K, y = `Pr[X|K]`), color = marker) +
    geom_line(aes(linetype = marker, color = marker)) +
    geom_point(aes(shape = marker, color = marker)) +
    theme_classic() +
    theme(strip.text.x = element_text(face = "italic"),
          strip.background = element_rect(color = "white"),
          panel.grid.major.x = element_line(
            size = .1, color = "grey", linetype = "dotted")) +
    scale_linetype_manual(values = ltym) +
    scale_color_manual(values = colo) +
    scale_shape_manual(values = pchm) +
    geom_errorbar(
      aes(x = K, ymin = `Pr[X|K]` - `sd Pr[X|K]`,
          ymax = `Pr[X|K]` + `sd Pr[X|K]`, color = marker)) +
    facet_wrap(marker ~ species, scales = "free_y")
  ggsave(paste0("data/final/bestk_metrics_",
                substitute(colo), ".pdf"), width = sizep, height = sizep / aspp)
}
ggplot_bestk(colo = colm)
ggplot_bestk(colo = colmbw)
