library(ggplot2)
#check convergence for all runs
#dataset#k
runs <- system("find data/final/ -name *stlog | grep -v usat", intern = T)#because usat failed
#runs <- system("find data/final/ -name *stlog", intern = T)
  p <- runs %>%
    #plyr::aaply(1, function(x){
    sapply(function(x){
    dataset <- gsub("^.*run_(.*)/.*$","\\1",x) %>% as.character()
    k <- gsub("^.*K([1-9]).*$", "\\1", x)
    lk <- system(paste0("cat ", x, " | grep [0,5][0]:"), intern = T) %>%
      .[-grep("--", .)] %>%
      gsub("^.*(-\\d+).*$", "\\1", .) %>%
      c(dataset, k, .)
  }) %>% t %>% as.data.frame

  #change cols to numeric
  for (i in 2:ncol(p)){
    p[, i] <- as.numeric(as.character(p[, i]))
  }; rm(i)

  pcol <- nlevels(p$V1) %>% sqrt %>% floor()
  it <- ncol(p) - 2

for (datset in 1:nlevels(p$V1)){
  levelp <- levels(p$V1)[datset]
  pd <- filter(p, V1 == levelp)
  pdf(paste0("data/final/convergence", levelp, ".pdf"))
  par(mfrow = c(3, 3))
  for (k in 2:max(pd$V2)){
    pdk <- filter(pd, V2 == k)
    y_lim <- range(pdk[, -c(1, 2)])
    plot(1, type = "n", xlim = c(1, it), ylim = y_lim,
         main = paste0("K = ", k), ylab = "Ln Like",
         xlab = "iterations x UPDATEFREQ (50) ")
    for (lin in 1:nrow(pdk)){
      lines(x = c(1:it),  y = pdk[lin, -c(1, 2)])
    }
  }
  dev.off()
}
#check for convergence
#(c) Jeff Hanson
#issue: read from parameters
#issue: be careful because for usatellites the run length is different
start_itertion_number <- 5050 #starting iteration after burnin
thinning_number <- 50 #frequency of saving steps

result <- dplyr::mutate_if(p, is.factor, as.character) %>%
 plyr::dlply(c("V1", "V2"), function(x) {
  y <- t(as.matrix(x[, c(-1, -2), drop = FALSE]))
  y <- as.numeric(y)
  y <- as.list(as.data.frame(matrix(y, nrow = ncol(x) - 2)))
  l <- lapply(y, function(z) {
    coda::mcmc(z, start = start_itertion_number,
               thin = thinning_number)
  })
  gr <- coda::gelman.diag(coda::mcmc.list(l), autoburnin = FALSE)
  data.frame(id = as.character(x$V1[1]), k = as.numeric(as.character(x$V2[1])),
             rhat = gr$psrf[[1]], upper_ci =  gr$psrf[[2]])
}
) %>% plyr::rbind.fill()
#rhat values above 1.05 imply chains have not converged.

pdf("data/final/rhat.pdf", height = 4, width = 6)
ggplot(data = result, aes(x = k, y = rhat, color = id)) + geom_line()
dev.off()
write.csv(result, file = "data/final/rhat.csv")
saveRDS(result, "data/intermediate/rhat.rds")
#run clumpak

#matrix correlation between:
# samples from same locality
# subsampled datasets:
#y = variance between runs
#x = k
#run best k
  #output table with bestk for all datasets
#structure maps for all datasets for best K, or alternatively if
#no best k, the choose K == 2
