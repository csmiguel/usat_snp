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

#load kfinder results
kfinder_ls <- readRDS("data/intermediate/kfinder_ls.rds")
source("code/functions/plot_kfinder.r")

#replace names
names(kfinder_ls) <-
  plyr::mapvalues(
    x = names(kfinder_ls),
    from = c("run_1000", "run_10000", "run_200", "run_20000", "run_3000",
      "run_500", "run_5000", "run_dart_hyla", "run_dart_pelo",
      "run_usat_hyla", "run_usat_pelo"),
    to = c("subset 1000", "subset 10000", "subset 200", "subset 20000",
           "subset 3000", "subset 500", "subset 5000", "SNPs H. molleri",
           "SNPs P. cultripes", "microsatellites H. molleri",
           "microsatellites P. cultripes"))

#plot parsimony

pdf("data/final/kfinder_subsampling.pdf", height = 9, width = 13)
  par(mfrow = c(3, 3))
  for (i in names(kfinder_ls[[1]])[-1]){
    plot_k(klist = kfinder_ls, index = i,
      runs_pattern = "[0-9]", mode = "subsampling")
    plot_k(klist = kfinder_ls, index = i, runs_pattern = "cult")
    plot_k(klist = kfinder_ls, index = i, runs_pattern = "molle")
  }
dev.off()
