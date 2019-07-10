nboot = 1000 #number of bootstrap
thresh = .7 #treshold for coloring nodes
#subsample sizes for pelobates dartseq data
s <- c(20000, 10000, 5000, 3000, 1000, 500, 200)
#times the markers are subsampled from the original matrix of genotypes
n_sampling_replicas <- 5
#maf thresholds
maf <- seq(0.04, 0.1, 0.02)
