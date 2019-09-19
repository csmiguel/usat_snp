###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# July 2019
###.............................................................................
#GOAL: get lambda and check convergence
#DESCRIPTION:
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

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

#plot convergence
for (datset in 1:nlevels(p$dataset)){ #for each dataset
  levelp <- levels(p$dataset)[datset] #name of dataset
  pd <- filter(p, dataset == levelp) #filter dataset
  pdf(paste0("data/final/convergence", levelp, ".pdf"))
  par(mfrow = n2mfrow(length(unique(pd$k))))
  for (K in min(pd$k):max(pd$k)){
    pdk <- filter(pd, K == k)
    y_lim <- range(pdk[, -c(1, 2)])
    plot(1, type = "n", xlim = c(1, it), ylim = y_lim,
         main = paste0("K = ", K), ylab = "Ln Like",
         xlab = "iterations x UPDATEFREQ (50) ")
    for (lin in 1:nrow(pdk)){
      lines(x = c(1:it),  y = pdk[lin, -c(1, 2)])
    }
  }
  dev.off()
}
#save lambda
saveRDS(lambda, "data/intermediate/lambda.rds")
