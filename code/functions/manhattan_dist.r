#calculate manhattan distance per locus
#ONLY use for microsatellite or small matrices of SNP data.
#For SNP data it is more efficient to use directly: stats::dist(l@tab, "manhattan")
#Why I developped this function?
# because for usat data there is no package to calculate Manhattan distances
# I dealt manually with missing data as in stats::dist
# It has been tested only with diploids
# The core of this function calculates Manhattan distance for each loci and then
# corrects for missing data.

per.loc <- function(l) {
  # I will divide by 2 because it counts alleles differences twice (1 per allele column)
  stats::dist(l@tab, "manhattan") / 2
}

genind_manhattan <- function(x){
  library(dartR)
  library(dplyr)
  assertthat::assert_that(adegenet::is.genind(x))
  assertthat::assert_that(all(ploidy(x) == 2))
  #distance per locus
  dists <- lapply(adegenet::seploc(x), per.loc)
  #data frame with distances
  hh <- plyr::ldply(dists, as.numeric) %>% .[, -1]
  #total distances per pairwise comparison
  hhd <- hh %>% apply(2, sum, na.rm = T)
  #total number of non-NA comparisons
  hhna <- apply(!is.na(hh), 2, sum)
  #no. of differences * (total loci / non-na loci)
  v <- hhd * (adegenet::nLoc(x) / hhna)
  #format data as distance matrix
  dist.mat <- matrix(NA, ncol = adegenet::nInd(x), nrow = adegenet::nInd(x))
  dist.mat[which(lower.tri(dist.mat) == TRUE)] <- v
  dist.mat <- as.dist(dist.mat)
  attr(dist.mat, "Labels") <- adegenet::indNames(x)
  attr(dist.mat, "method") <- "Manhattan"
  attr(dist.mat, "call") <- "stats::dist per locus"
  return(dist.mat)
}
