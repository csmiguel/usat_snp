###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: filter datasets of SNPs and usats to have the same samples
#DESCRIPTION:
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

library(dartR)
library(dplyr)
#which samples are shared between datasets:
gen <- readRDS("data/intermediate/gen_consolidated_filtered.rds")

# remove individuals from genind obj
# vector with individuals to keep
samples_overlapped_dartusat <-
  indNames(gen$dart_pelo)[(adegenet::indNames(gen$dart_pelo) %in% adegenet::indNames(gen$usat_pelo))]
# new DART genind object
overlapped_dart_pelo <-
  gen$dart_pelo[row.names(gen$dart_pelo@tab) %in% samples_overlapped_dartusat]
# new usat genind object
overlapped_usat_pelo <-
  gen$usat_pelo[row.names(gen$usat_pelo@tab) %in% samples_overlapped_dartusat]

# check number of Individuals is
assertthat::assert_that(nInd(overlapped_usat_pelo) == nInd(overlapped_dart_pelo))

gen <- c(overlapped_dart_pelo, overlapped_usat_pelo)
names(gen) <- c("dart_pelo", "usat_pelo")

saveRDS(gen, "data/intermediate/shared_gen_consolidated_filtered.rds")
