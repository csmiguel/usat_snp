ind_miss_thresh <- 0.65 #missingness allowed
repavg_threshold <- 0.95
locus_callrate_threshold <- 0.8 #allows 1-locus_callrate_threshold missingness
sec_method <- "best"
min_maf <- 0.02 #eg: 95samples*2ploidy*0.8callrate*0.02maf = 3 alleles
hwe_alpha <- 0.05
cov_th <- 2 #treshold for mean diff. in coverage between ref and alt alleles
balance_tresholds <- c(0.15, 0.85)
max_coverage <- 3.5 # deviations from the median
