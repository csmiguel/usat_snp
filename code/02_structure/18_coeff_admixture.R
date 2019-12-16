###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: compute coeffecient of admixture. This script computes for the Pelobates
# dataset alone a coefficient of admixture. This coeffecient is computed for
# each individual as the sum of square ancestry memberships for each K corrected
# for the maximum possible value that it could reach (if memberships across all
# K's were the same. It substracts the values between marker types. So, positive
# values indicate Dart genotypes yield STRUCTURE ancestries which are less
# admixed while negative results indicate the opposite. The absoluste numeric
# value can be interpreted as the proportion of how much more admixed is one
# dataset respect to the other, being 0 equal and 1 the maximum possible.
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(dplyr)
library(ggplot2)

qclumpak <- readRDS("data/intermediate/clumpak_major.rds")
source("code/functions/coeff_admixture.r")

#select Pelobates data
qclumpak <- qclumpak[grep(pattern = c("usat_pelo|dart_pelo"), x = names(qclumpak))]

# pairwise computes admixture coefficient
qpairwise <-
  lapply(seq_along(qclumpak[[1]]), function(k){
    husat <- apply(qclumpak[["usat_pelo"]][[k]], 1, ca)
    hdart <- apply(qclumpak[["dart_pelo"]][[k]], 1, ca)
    h <- hdart - husat
  }) %>% reshape2::melt() %>%
  dplyr::filter(L1 > 1)
names(qpairwise) <- c("ca", "k")

saveRDS(qpairwise, "data/intermediate/coeff_admixture.rds")
