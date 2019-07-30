###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: plot black lines convergence of replicate runs across K for each dataset
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(dplyr)

conv <- readRDS("data/intermediate/prob_convergence.rds")
source("code/functions/read_lk_str.r")
source("code/parameters/structure.r")

#Plot convergence
plot_convergence(method = "with_burnin", p = conv$usat, burnin = burnin_usat)
plot_convergence(method = "with_burnin", p = conv$dart, burnin = burnin_snp)
