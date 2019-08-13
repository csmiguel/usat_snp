###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: check convergence for structure runs
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

library(dplyr)
source("code/functions/read_lk_str.r")
source("code/parameters/structure.r")

#create vector with path for structure log files
runs <- run("all")
runs <- runs[grep(pattern = "shared", runs)]
#read all log files into a dataframe with:#rowname, path of file#dataset,
#dataset#k, k# columns with names being the interation and value being the lnlk.
pusat <- runs %>% grep(pattern = "usat", value = T) %>% read_lk()
pdart <- runs %>% .[-grep(pattern = "usat", .)] %>% read_lk()

#check for convergence
#(c) Jeff Hanson
convergence_usat <- rhat(p = pusat, burnin = burnin_usat)
convergence_dart <- rhat(p = pdart, burnin = burnin_snp)

#write results
write.csv(convergence_usat, file = "data/final/shared_convergence_usat.csv")
write.csv(convergence_dart, file = "data/final/shared_convergence_dart.csv")

saveRDS(list(usat = convergence_usat, dart = convergence_dart),
  "data/intermediate/shared_rhat.rds")
saveRDS(list(usat = pusat, dart = pdart),
  "data/intermediate/shared_prob_convergence.rds")
