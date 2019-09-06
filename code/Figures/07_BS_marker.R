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

plot_data <- readRDS("data/intermediate/bs_treemetricsPlot.rds")
pred_data <- readRDS("data/intermediate/model_marker_pred_data.rds")

source("code/parameters/plotting_par.r")
source("code/functions/plot_bs.r")

#labels distances
labels_distance <- c(dnp = "Dnp", dnr = "Dnr")
labels_species <- c(hyla = "H. molleri", pelo = "P. cultripes")

#produce plots
plot_bs(mode = "bw")
plot_bs(mode = "color")
