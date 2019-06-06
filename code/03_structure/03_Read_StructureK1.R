library(ggplot2)
#check convergence for all runs
#dataset#k
source("code/functions/read_lk_str.r")
runs <- run("k1")
#read all log files into a dataframe with:#rowname, path of file#dataset,
#dataset#k, k# columns with names being the interation and value being the lnlk.
p <- read_lk(runs = runs)
# get lambda
lambda <- get_lambda()

  pcol <- nlevels(p$V1) %>% sqrt %>% floor()
  it <- ncol(p) - 2

for (datset in 1:nlevels(p$V1)){#for each dataset
  levelp <- levels(p$V1)[datset] #name of dataset
  pd <- filter(p, V1 == levelp)#filter dataset
  pdf(paste0("data/final/convergence", levelp, ".pdf"))
  par(mfrow = n2mfrow(length(unique(pd$V2))))
  for (k in min(pd$V2):max(pd$V2)){
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
saveRDS(lambda, "data/intermediate/lambda.rds")
