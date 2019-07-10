# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: Plot results from subsampling SNPs from Pelobates
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(dplyr)
library(magrittr)

gen <- readRDS("data/intermediate/gen_consolidated_filtered.rds")
nloc_pelo <- gen$dart_pelo %>% adegenet::nLoc()

bstree <- readRDS("data/intermediate/s_bs_treemetrics.rds") %>%
  {c(list(
    dart_pelo = readRDS("data/intermediate/bs_treemetrics.rds")$dart_pelo), .)}

bt <- readRDS("data/intermediate/bs_treemetrics.rds")
bt %<>% .[!grepl("dart_pelo", names(bt))]

names(bstree)[names(bstree) == "dart_pelo"] <- nloc_pelo

# plot parameters
x_lab <- c("Distance to parent node", "Distance from root")
max_dnr <- sapply(bstree, function(x) x$dnr) %>% unlist %>% max
max_dnp <- sapply(bstree, function(x) x$dnp) %>% unlist %>% max
x_lim <- matrix(c(0, max_dnp, 0, max_dnr), ncol = 2)

#1. Point plots

pdf(file = paste0("data/final/", "subsampling_bs", ".pdf"),
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
    legend("bottomright", legend = paste(names(bstree), "loci"),
           pch = 20, col = grey.colors(length(bstree)), bty = "n")
    }
  dev.off()

#2. Boxplot
#labs for plot
nl <- sapply(gen, adegenet::nLoc) %>% paste("loci")

#dataframe where col1 is bootstrap support and col2 is dataset
s_bt_m <-
  c(bt, bstree) %>% lapply(function(x) dplyr::select(x, bs_tm)) %>%
  reshape2::melt() %>%
  dplyr::select(-variable) %>%
  dplyr::rename(dataset = 2) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(dataset = as.factor(dataset))
#reorder levels
s_bt_m$dataset <-
  factor(s_bt_m$dataset,
         levels = c("200", "500", "1000", "3000", "5000", "10000", "20000",
                    "33144", "usat_pelo", "dart_hyla", "usat_hyla"))

lab <- paste(levels(s_bt_m$dataset),
  nl[match(levels(s_bt_m$dataset), names(gen))]) %>%
  gsub(pattern = " NA", replacement = "") %>%
  gsub(pattern = "usat_hyla", replacement = "microsatellites \nH. molleri,") %>%
  gsub(pattern = "usat_pelo",
    replacement = "microsatellites \nP. cultripes,") %>%
  gsub(pattern = "dart_hyla", replacement = "SNPs \nH. molleri,")
colbp <- rep("grey", nlevels(s_bt_m$dataset))
colbp[grep("usat", levels(s_bt_m$dataset))] <- grey(0.3)

#boxplot
pdf(file = "data/final/bs_boxplot.pdf", height = 5.5, width = 8)
par(mar = c(7, par()$mar[2:4]))
boxplot(s_bt_m$value~s_bt_m$dataset, xlab = "Number of loci",
        ylab = "Bootstrap support",
        col = colbp,
        pch = 20, outcol = grey(0.7), xaxt = "n", yaxt = "n")
axis(1, at = seq_along(levels(s_bt_m$dataset)), labels = lab,
  las = 2, cex.axis = 0.8)
axis(2, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), las = 1, cex.axis = 0.8)
abline(v = length(bstree) - 0.5, lty = 2, lwd = 2)
dev.off()
