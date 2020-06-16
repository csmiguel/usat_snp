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
library(dplyr)

qpairwise <-
  qpairwise %>%
  dplyr::mutate(k = as.factor(k)) %>%
  dplyr::as_tibble()

ggplot(qpairwise, aes(y = k, x = ca, fill = ..x..)) +
  geom_density_ridges_gradient() +
  scale_fill_gradient2(low = "red",
                          mid = "white",
                          high = "blue",
                       name = "admixture") +
  xlab("Coefficient of admixture") +
  ylab("K") +
  theme_ridges() +
  scale_y_discrete(limits = c(rev(levels(qpairwise$k)))
                  ) +
  scale_x_continuous(limits = c(-1, 1),
                     expand = c(0, 0),
                     position = "top") +
  geom_vline(xintercept = 0, linetype = 2, color = "blue", size = 0.5)

ggsave("data/final/coeff_admixture.pdf", height = 4, width = 6)
