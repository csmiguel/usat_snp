#create data frame with combined BS for subsampled trees and original ones
s.btx <- lapply(s.bt, function(x) x[-1])
s_bt_m <- reshape::melt(s.btx)
s_bt_m[, 1] <- s_bt_m[, 1] / (nboot * n_sampling_replicas)
assertthat::assert_that(max(s_bt_m[, 1]) == 1)
s_bt_m[, 2] <- as.numeric(s_bt_m[, 2]) %>% as.factor()
s_bt_m <- rbind(s_bt_m, reshape::melt(bt))
names(s_bt_m) <- c("BS", "dataset")
s_bt_m$BS <- round(s_bt_m$BS, 2)

#quantiles for BS
g <- by(data = s_bt_m$BS, INDICES = s_bt_m$dataset, quantile)
gg <- 1:length(g) %>% sapply(function(x) g[[x]] %>% round(2))
colnames(gg) <- names(g); rm(g)

#labs for plot
nl <- sapply(gen, adegenet::nLoc)
lab <- c(levels(s_bt_m$dataset)[1:length(s)],
  paste0(names(gen), rep("\n (", 4), nl, rep(")", 4)))

#boxplot
pdf(file = "data/intermediate/bs_boxplot.pdf", height = 5, width = 8)
  boxplot(s_bt_m$BS~s_bt_m$dataset, xlab = "Number of loci",
      ylab = "Bootstrap support",
      main = "Change of overall boostrap support accross number of loci",
      col = c(rep("grey", length(s)), rep(grey(0.3), 2), rep(grey(0.8), 2)),
      pch = 20, outcol = grey(0.7), xaxt = "n", yaxt = "n")
      axis(1, at = 1:(length(s) + 4), labels = lab, las = 2)
      axis(2, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), las = 1)
      abline(v = length(s) + 0.5, lty = 2, lwd = 2)
dev.off()


#save trees (first tree is the real tree)
saveRDS(s.tr, "data/intermediate/s_trees.rds")
saveRDS(gg, "data/intermediate/s_quantiles.rds")
