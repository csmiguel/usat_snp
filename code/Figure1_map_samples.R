###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: Plot map wih samples
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
#plotting function
source("code/functions/map_distribution.r")
#load consolidated genotypes
gen <- readRDS("data/intermediate/gen_consolidated_filtered.rds")

plot_samples(meta = gen,
  name_plot = "data/final/Figure1_map_samples.pdf", res = 50000)
