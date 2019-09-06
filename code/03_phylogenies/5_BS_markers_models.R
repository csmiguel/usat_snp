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
# Initialization
## set paths
rm(list = ls())
data_path <- "~/Downloads/bs_treemetrics.rds"
source("~/Downloads/plotting_par.r")

## settings
standard_error_multiplier <- 1.96

## load packages
library(ggplot2)
library(dplyr)
library(tibble)
library(broom)
library(tidyr)

# Preliminary processing
## load data
raw_data <- readRDS(data_path)

## add columns to data
marker <- sapply(strsplit(names(raw_data), "_", fixed = TRUE), `[[`, 1)
species <- sapply(strsplit(names(raw_data), "_", fixed = TRUE), `[[`, 2)
data <- plyr::ldply(seq_along(raw_data), function(i) {
  x <- raw_data[[i]]
  x$marker <- marker[i]
  x$species <- species[i]
  x
}) %>%
as_tibble()

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
  suppressWarnings(tidy(anova(x$full, x$null, test = "Chisq")))
})
dnr_test <- lapply(dnr_models, function(x) {
  suppressWarnings(tidy(anova(x$full, x$null, test = "Chisq")))
})

# Exports
## format raw data for plotting
plot_data <-
  data %>%
  gather(dist, value, -bs_tm, -marker, -species)

## format model predictions for plotting
### make predictions for dnp
pred_dnp_data <-
  plyr::ldply(dnp_models, function(x) {
    d <- expand.grid(dnp = seq(0, max(data$dnp), length.out = 1000),
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
    d <- expand.grid(dnr = seq(0, max(data$dnr), length.out = 1000),
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
pred_data <- bind_rows(pred_dnp_data, pred_dnr_data)

## make plot
labels_distance <- c(dnp = "Dnp", dnr = "Dnr")
p <-
  ggplot() +
  geom_ribbon(aes(x = value, ymin = bs_tm_lower, ymax = bs_tm_upper,
                  fill = marker), size = 1e-5, alpha = 0.1, data = pred_data) +
  geom_point(aes(x = value, y = bs_tm, colour = marker, shape = marker),
             data = plot_data) +
  geom_line(aes(x = value, y = bs_tm, color = marker), data = pred_data) +
  facet_grid(species ~ dist, scales = "free_x",
             labeller = labeller(dist = labels_distance)) +
  labs(y = "Bootstrap support", x = "") +
  theme(strip.text.y = element_text(size = 10, face = "italic"),
        strip.text.x = element_text(size = 10)) +
  scale_color_manual(values = c("black", "red")) +
  scale_fill_manual(values = c("black", "red")) +
  scale_shape_manual(values = c(16, 16)) +
  theme(legend.position = c(0.35, 0.3),
        legend.background = element_rect(size = 0.2, linetype = "solid",
                                         colour = "black"),
        legend.text = element_text(colour = "black", size = 10),
        legend.title =  element_blank())

## render plot
print(p)
