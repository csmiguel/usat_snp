###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: plot
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(ggplot2)
best_files <- readRDS("data/intermediate/best_files_kfinder.rds")
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
