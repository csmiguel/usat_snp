#plot results from bootstrap analysis
bstree <- readRDS("data/intermediate/bs_treemetrics.rds")
# plot parameters
x_lab <- c("Distance to parent node", "Distance from root")
max_dnr <- sapply(bstree, function(x) x$dnr) %>% unlist %>% max
max_dnp <- sapply(bstree, function(x) x$dnp) %>% unlist %>% max
x_lim <- matrix(c(0, max_dnp, 0, max_dnr), ncol = 2)

#produce plots
for (i in 1:2){
  pdf(file = paste0("data/final/",
  gsub(" ", "_", x_lab[i]), ".pdf"),
   height = 8, width = 10)
  par(mfrow = c(2, 2))
    for (p in 1:length(bstree)){
      hh <- bstree[[p]]
      plot(as.numeric(as.matrix(hh[, i])),
          hh$bs_tm, pch = 20,
          xlab = x_lab[i],
          ylab = "Bootstrap support",
          main = names(bstree)[p],
          xlim = x_lim[, i],
          ylim = c(0, 1))
      }
  dev.off()
}
