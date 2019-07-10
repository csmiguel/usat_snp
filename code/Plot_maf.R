# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: Plot results from maf filtering in SNPs from Pelobates
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

library(dplyr)
library(magrittr)

bstree <- readRDS("data/intermediate/maf_bs_treemetrics.rds") %>% rev

# plot parameters
x_lab <- c("Distance to parent node", "Distance from root")
max_dnr <- sapply(bstree, function(x) x$dnr) %>% unlist %>% max
max_dnp <- sapply(bstree, function(x) x$dnp) %>% unlist %>% max
x_lim <- matrix(c(0, max_dnp, 0, max_dnr), ncol = 2)

#1. Point plots

pdf(file = paste0("data/final/", "maf_bs", ".pdf"),
    height = 9, width = 6.5)
par(mfrow = c(2, 1), las = 1)
for (i in 1:2){
  plot(1, pch = 20,
       xlab = x_lab[i],
       ylab = "Bootstrap support",
       xlim = x_lim[, i],
       ylim = c(0, 1), type = "n")
  for (p in seq_along(bstree)){
    hh <- bstree[[p]]
    points(as.numeric(as.matrix(hh[, i])),
           hh$bs_tm, pch = 20, col = grey.colors(length(bstree))[p])
  }
  legend("bottomright", legend = names(bstree),
         pch = 20, col = grey.colors(length(bstree)), bty = "n")
}
dev.off()

#2. Boxplot

#dataframe where col1 is bootstrap support and col2 is dataset
s_bt_m <-
  bstree %>% lapply(function(x) dplyr::select(x, bs_tm)) %>%
  reshape2::melt() %>%
  dplyr::select(-variable) %>%
  dplyr::rename(dataset = 2) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(dataset = as.factor(dataset))

#boxplot
pdf(file = "data/final/maf_bs_boxplot.pdf", height = 5.5, width = 8)
par(las = 1)
boxplot(s_bt_m$value~s_bt_m$dataset, xlab = "MAF threshold",
        ylab = "Bootstrap support",
        col = "grey",
        pch = 20, outcol = grey(0.7))

dev.off()

#3. Point plot linear model
h <- s_bt_m
# model
m1 <- lm(value ~ as.numeric(as.character(dataset)), data = h)
summary(m1)
pr <- predict(m1,
  newdata = data.frame(dataset = seq(0, 0.1, length.out = 1000)))

pdf(file = "data/final/maf_model.pdf", height = 4.5, width = 6)
plot(seq(0, 0.1, length.out = 1000),
     pr,
     ylab = "Bootstrap support",
     xlab = "MAF", las = 1, type = "l", ylim = c(0, 1))
points(jitter(as.numeric(as.character(h$dataset))),
       h$value,
       col = adjustcolor("blue", alpha.f = 0.2), pch = 20)
dev.off()
