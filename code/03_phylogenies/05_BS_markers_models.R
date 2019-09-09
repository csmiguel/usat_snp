###.............................................................................
# (c) Jeffrey Hanson developped most of this script
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: fit models for relation between BS and markers (for each sp and metric)
#DESCRIPTION:
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

library(dplyr)
library(tibble)
library(broom)
library(tidyr)

raw_data <- readRDS("data/intermediate/bs_treemetrics.rds")

## settings
standard_error_multiplier <- 1.96

# Preliminary processing

## prepare data
marker <- sapply(strsplit(names(raw_data), "_", fixed = TRUE), `[[`, 1)
species <- sapply(strsplit(names(raw_data), "_", fixed = TRUE), `[[`, 2)
data <- plyr::ldply(seq_along(raw_data), function(i) {
  x <- raw_data[[i]]
  x$marker <- marker[i]
  x$species <- species[i]
  x
}) %>%
as_tibble()

# 1. Compare BS averages between makers for each species

bs_models <- plyr::dlply(data, c("species"), function(x) {
  suppressWarnings({
    glm(bs_tm ~ marker, data = x, family = "binomial") %>% summary()
  })
})

sink("data/final/BS_models_means_marker.txt")
bs_models
sink()

# 2. Compare change of Dnp and Dnr between markers for each species

# Main processing
## fit models
dnp_models <- plyr::dlply(data, c("species"), function(x) {
  suppressWarnings({
    full <- glm(bs_tm ~ dnp * marker, data = x, family = "binomial")
    null <- update(full, . ~ . - dnp:marker)
  })
  list(full = full, null = null)
})
dnr_models <- plyr::dlply(data, c("species"), function(x) {
  suppressWarnings({
    full <- glm(bs_tm ~ dnr * marker, data = x, family = "binomial")
    null <- update(full, . ~ . - dnr:marker)
  })
  list(full = full, null = null)
})

## test for difference between marker types
dnp_test <- lapply(dnp_models, function(x) {
  suppressWarnings(broom::tidy(anova(x$full, x$null, test = "Chisq")))
})
dnr_test <- lapply(dnr_models, function(x) {
  suppressWarnings(broom::tidy(anova(x$full, x$null, test = "Chisq")))
})

# Exports
## format raw data for plotting
plot_data <-
  data %>%
  tidyr::gather(dist, value, -bs_tm, -marker, -species)

## format model predictions for plotting
### make predictions for dnp
pred_dnp_data <-
  plyr::ldply(dnp_models, function(x) {
    d <- expand.grid(dnp = seq(0, max(x$full$data$dnp),
                     length.out = 1000),
                     marker = unique(data$marker))
    p <- predict(x$full, newdata = d, type = "link", se.fit = TRUE)
    d$bs_tm <- plogis(p$fit)
    d$bs_tm_lower <- plogis(p$fit - (standard_error_multiplier * p$se.fit))
    d$bs_tm_upper <- plogis(p$fit + (standard_error_multiplier * p$se.fit))
    d
  }) %>%
  rename(value = dnp) %>%
  mutate(dist = "dnp") %>%
  as_tibble()

### make predictions for dnr
pred_dnr_data <-
  plyr::ldply(dnr_models, function(x) {
    d <- expand.grid(dnr = seq(0, max(x$full$data$dnr), length.out = 1000),
                     marker = unique(data$marker))
    p <- predict(x$full, newdata = d, type = "link", se.fit = TRUE)
    d$bs_tm <- plogis(p$fit)
    d$bs_tm_lower <- plogis(p$fit - (standard_error_multiplier * p$se.fit))
    d$bs_tm_upper <- plogis(p$fit + (standard_error_multiplier * p$se.fit))
    d
  }) %>%
  rename(value = dnr) %>%
  mutate(dist = "dnr") %>%
  as_tibble()

### combine data into a single table
pred_data <- dplyr::bind_rows(pred_dnp_data, pred_dnr_data)

#export objects
saveRDS(pred_data, "data/intermediate/model_marker_pred_data.rds")
saveRDS(plot_data, "data/intermediate/bs_treemetricsPlot.rds")

#sink results models
sink("data/final/BS_models_dnpr_marker.txt")

cat("Results from models:\n1.1. Dnr models")
lapply(dnr_models, function(x) lapply(x, summary))
cat("\n\n1.2. Dnr test")
dnr_test
cat("\n\n2.1. Dnp models")
lapply(dnp_models, function(x) lapply(x, summary))
cat("\n\n2.2. Dnp test")
dnp_test
sink()
