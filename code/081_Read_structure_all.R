library(ggplot2)
#check convergence for all runs
#dataset#k
source("code/functions/read_lk_str.r")
source("code/parameters/structure.r")

#create vector with path for structure log files
runs <- run("all")
#read all log files into a dataframe with:#rowname, path of file#dataset,
#dataset#k, k# columns with names being the interation and value being the lnlk.
p <- read_lk()
  #change cols to numeric
  for (i in 2:ncol(p)){
    p[, i] <- as.numeric(as.character(p[, i]))
  }; rm(i)

#vector with iterations
it <- as.numeric(names(p)[3:length(names(p))])
#columns index which are burnin
burnin_cols <- c(1: which(names(p) == burnin_snp))

plot_convervenge(method = "with_burnin")
plot_convervenge(method = "no_burnin")


#check for convergence
#(c) Jeff Hanson
#issue: be careful because for usatellites the run length is different

thinning_number <- it[2] - it[1] #frequency of saving steps
#starting iteration after burnin
start_itertion_number <- burnin_snp + thinning_number

result <- p %>% select(-burnin_cols[-c(1, 2)]) %>%
  dplyr::mutate_if(is.factor, as.character) %>%
  plyr::dlply(c("dataset", "k"), function(x) {
  y <- t(as.matrix(x[, -c(1, 2), drop = FALSE]))
  y <- as.numeric(y)
  y <- as.list(as.data.frame(matrix(y, nrow = ncol(x) - 2)))
  l <- lapply(y, function(z) {
    coda::mcmc(z, start = start_itertion_number,
               thin = thinning_number)
  })
  gr <- coda::gelman.diag(coda::mcmc.list(l), autoburnin = FALSE)
  data.frame(id = as.character(x$dataset[1]),
             k = as.numeric(as.character(x$k[1])),
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
