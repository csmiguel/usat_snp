###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: compute standarized multilocus heterozygosity
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(dplyr)
library(ggplot2)

#load sMLH
het <- readRDS("data/intermediate/sMLH.rds") %>%
# select data for P. cultripes only
  {.[grep("pelo", names(.))]}

assertthat::assert_that(
  all(names(het[["usat_pelo"]]) == names(het[["dart_pelo"]])))
#pearson correlation
p1.test <- cor.test(het[["usat_pelo"]], het[["dart_pelo"]], method = "pearson")

#linear model
m1 <- lm(het[["dart_pelo"]] ~ het[["usat_pelo"]])

sink("data/final/correlation_sMLH.txt")
p1.test
sink()
#save residuals from linear model
saveRDS(m1$residuals, "data/intermediate/residuals-sMLH.rds")
