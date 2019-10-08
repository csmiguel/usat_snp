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
#select Pelobates data
qclumpak <- qclumpak[grep(pattern = c("usat_pelo|dart_pelo"), x = names(qclumpak))]

# computes admixture coefficient
qpairwise <-
  lapply(seq_along(qclumpak[[1]]), function(k){
    husat <- apply(qclumpak[["usat_pelo"]][[k]], 1, function(y) sum(y^2))
    hdart <- apply(qclumpak[["dart_pelo"]][[k]], 1, function(y) sum(y^2))
    h <- hdart - husat
  }) %>% reshape2::melt() %>%
  dplyr::filter(L1 > 1)

qpairwise$value_corrected <-
  apply(qpairwise, 1, function(x) x[1] / (1 - sum(rep(1 / x[2], x[2]) ^ 2)))

saveRDS(qpairwise, "data/intermediate/coeff_admixture.rds")
