###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: calculate best K
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
#  REQUIRED FILES:
#   Description:
#   Inpath:
#  OUTPUT:
#    Description:
#    Outpath:
#  DEPENDENCIES:
###.............................................................................
library(dplyr)
library(magrittr)
library(dplyr)
#1. create parameters file
# vector with folder names for each run
kfinder <- dir("data/final", pattern = "str.K$", full.names = T, recursive = T)
#create list of data frames with K statistics
kfinder_ls <-
  kfinder %>%
  lapply(function(x){
    readLines(x) %>%
      grep(pattern = "^\\s+[0-9]", value = T) %>%
      textConnection() %>%
      read.table() %>%
      setNames(c("K", "Pritchard", "Evanno", "Parsimony"))
  })
names(kfinder_ls) <- kfinder %>%
  stringr::str_extract("run.*/") %>%
  gsub(pattern = "/", replacement = "")

#plot parsimony
plot_k <- function(klist = NULL, index = NULL, runs_pattern = NULL, mode = 1){
  #runs_pattern, is the pattern to grep from names
  #index is Pritchard, Evanno or Parsimony
  #list with data frames of k statistics
  to_plot <- klist[grep(runs_pattern, names(klist))]
  if (mode == "subsampling"){
    to_plot <-
      to_plot[names(to_plot) %>%
                stringr::str_extract("[0-9]+") %>% as.numeric() %>% order]
  }
  y_lim <- to_plot %>% sapply(function(x) x[[index]]) %>%
    unlist %>% as.numeric() %>% range(na.rm = T)
  plot(1, type = "n",
       xlim = range(to_plot[[1]]$K),
       ylim = y_lim,
       ylab = index,
       xlab = "K")
  cols <- colorRampPalette(c("blue", "red"))
  for (i in seq_along(to_plot)){
    lines(to_plot[[i]]$K,
          as.numeric(as.character(to_plot[[i]][[index]])),
          lwd = 2,
          type = "o",
          col = cols(length(to_plot))[i]
          )
  }
  legend("bottomright", legend = names(to_plot),
         col = cols(length(to_plot)), pch = 1, lwd = 2,
         bg = adjustcolor("white", alpha.f = 0.5))
}
plot_k(klist = kfinder_ls, index = "Parsimony", runs_pattern = "pelo")
plot_k(klist = kfinder_ls, index = "Evanno", runs_pattern = "hyla")
plot_k(klist = kfinder_ls, index = "Parsimony", runs_pattern = "[0-9]", mode = "subsampling")




#best FILES
best_files <-
  kfinder %>%
  sapply(function(x){
    readLines(x) %>%
      grep(pattern = "Best output", value = T) %>%
      sub(pattern = "^.* /", replacement = "/", .)
  })

#read structure output
h <- pophelper::readQStructure(best_files[1], indlabfromfile = T)

#sort indlabels by populations
# read genotypes
gen <- readRDS("data/intermediate/gen_consolidated_filtered.rds")
pelo <- gen$dart_pelo$other$metadata %>%
          dplyr::arrange(locality)

#sort ancestry table according to order of populations
h$str_K6_rep3_f %<>%
  .[match(rownames(.), pelo$sample_id), ]
#plot structure output
p <-
  pophelper::plotQ(h,
                   exportplot = F,returnplot = T,
                   showindlab = T, indlabsize = 8,
                   height = 3, barsize = 1, useindlab = T,
                   showtitle = T, barbordercol = "white",
                   showyaxis = T, sortind = "label",
                   ordergrp=T, grpmean = F, grplabsize=2,
                   grplabangle = 60,linepos = 1,
                   grplabpos=0.9,splabsize = 8,
                   grplabheight = 3)

plot(p$plot[[1]])
