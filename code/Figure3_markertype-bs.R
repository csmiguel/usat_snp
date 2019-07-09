###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: Plot results from bootstrap analysis
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
#plot results from bootstrap analysis
bstree <- readRDS("data/intermediate/bs_treemetrics.rds")
# plot parameters
x_lab <- c("Distance to parent node", "Distance from root")
max_dnr <- sapply(bstree, function(x) x$dnr) %>% unlist %>% max
max_dnp <- sapply(bstree, function(x) x$dnp) %>% unlist %>% max
x_lim <- matrix(c(0, max_dnp, 0, max_dnr), ncol = 2)

#rename datasets
bstree
names(bstree) <-
  plyr::mapvalues(
    x = names(bstree),
    from = names(bstree),
    to = c("SNPs H. molleri", "SNPs P. cultripes",
    "microtellites H. molleri", "microtellites P. cultripes"))
#produce plots
for (i in 1:2){
  pdf(file = paste0("data/final/",
  gsub(" ", "_", x_lab[i]), ".pdf"),
   height = 7, width = 8.8)
  par(mfrow = c(2, 2), las = 1)
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
