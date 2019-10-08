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

qpairwise <- readRDS("data/intermediate/coeff_admixture.rds")
library(ggplot2)

ggplot(qpairwise) +
  geom_boxplot(aes(y = value_corrected, x = as.factor(L1)), width = 0.1) +
  theme_classic() +
  ylab("Coefficient of admixture") +
  xlab("K") +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
  geom_abline(slope = 0, intercept = 0, linetype = 2, color = "blue", size = 0.3)


  ggsave("data/final/coeff_admixture.pdf", height = 4, width = 6)
