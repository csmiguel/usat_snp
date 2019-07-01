#calculate difference in coverage. I will calculate the difference in coverage
# between alleles, units being number of times the most covered allele is
#different from the other. ie: 2.5 means that the allele with most converage
#has in average 2.5 times more coverage than the allele with less converage.
#This step was applied in 10.1111/1755-0998.12997 (Tasmanian devil)
#Thank Elspeth in acknowledgements for sharing part of the code.
diff_coverage <- function(genlight_object){
  snp <- genlight_object$other$loc.metrics
  xx <- 1 / (1 - ( (abs(snp[, "AvgCountRef"] - snp[, "AvgCountSnp"])) /
            pmax(snp[, "AvgCountRef"], snp[, "AvgCountSnp"])))
  assertthat::assert_that(length(xx) == adegenet::nLoc(genlight_object))
          }
filt_dartgenlight <- function(genlight, filt_vector){
  #this function removes the loci from a genlight object
  #together with the corresponding columns from $loc.metrics
  #the genlight object must be formatted with dartR::gl.read.dart
  #It only retains loci from filt_vector.
   if (length(filt_vector) != adegenet::nLoc(genlight))
    warning("length of filt_vector is different to no loci in genlight")
  assertthat::assert_that(nrow(other(genlight)$loc.metrics) ==
    adegenet::nLoc(genlight),
      msg = "number of loci do not match with nrows in loc.metrics")

  new_gen <- genlight[, filt_vector]
  adegenet::other(new_gen)$loc.metrics <-
    adegenet::other(new_gen)$loc.metrics[filt_vector, ]
  new_gen
}

all_balance_filter <- function(genlight, tresholds = c(0.2, 0.8),
include_plot = NULL, plot_name){
  #this function takes genlight formatted with darR::gl.read.dart function
  #and removes unbalanced alleles below and above the tresholds, respectively.
  #it can return plot in working dir. Based on publication 10.1111/mec.14792
  # O'Leary et al. 2018
  #usage:
  #all_balance_filter(gen, include_plot = T, plot_name = "Allele_balance_filt")
  assertthat::assert_that(class(genlight) == "genlight")
  snp <- genlight$other$loc.metrics
  ab <- snp$AvgCountSnp / (snp$AvgCountSnp + snp$AvgCountRef)
  removal_vector <- which(ab > tresholds[1] & ab < tresholds[2])
if (include_plot == T){
  assertthat::assert_that(is.character(plot_name))
  pdf(file.path(getwd(), paste0(plot_name, ".pdf")), width = 5, height = 4)
  hist(ab, xlab = "Allele balance", main = "",
       col = adjustcolor("orange", 0.5), breaks = 20, xlim = c(0, 1))
  abline(v = c(tresholds[1], 0.5, tresholds[2]), lty = 2,
         col = c("red", "green", "red"), lwd = 2)
  dev.off()
  }
  output_gen <- filt_dartgenlight(genlight, removal_vector)
  cat("Filter on allele balance removed",
    adegenet::nLoc(genlight) - adegenet::nLoc(output_gen),
   "loci from a total of", adegenet::nLoc(genlight),
   "loci. The filtered genlight object has", adegenet::nLoc(output_gen),
   "loci, using thresholds:", tresholds,
   if (include_plot == T) ". Plot saved at",
   file.path(getwd(), paste0(plot_name, ".pdf\n")))
     return(output_gen)
}

#funtion to filter genlight object formatted with dartR::gl.read.dart function
#It removes loci that fall outside median_dev times from the median.
#See figure 4 in 10.1111/mec.14792.
#We do not know exactly how DArTseq calculates AvgCountRef and
#AvgCountSnp columns. Indeed, we do not know if they correct for sample
#sequencing depth or ploidy. For that reason, we choose a relaxed filtering.
#in 10.1111/mec.14792, fig4, they show an example with 2 times the mode, but
#no lower limit is included.
total_coverage_filtering <- function(genlight, median_dev, include_plot = NULL,
  plot_name, lower_end = NULL){
  snp <- genlight$other$loc.metrics %>%
    mutate(comb = AvgCountRef + AvgCountSnp) %>%
    dplyr::select(comb) %>% .[, 1]
  if (lower_end == T){
    thresholds <- c(median(snp) * median_dev, median(snp) * 1 / median_dev)
    filt_vector <- which(snp < thresholds[1] & snp > thresholds[2])
  } else if (lower_end == F){
    thresholds <- c(median(snp) * median_dev)
    filt_vector <- which(snp < thresholds[1])
  }
  output_gen <- filt_dartgenlight(genlight, filt_vector)
  #plot
  if (include_plot == T){
  colb <- adjustcolor("orange", alpha.f = 0.5)
  bz <- 3 # breaks
  pdf(file.path(getwd(), paste0(plot_name, ".pdf")), width = 5, height = 4)
  hist(snp, main = "",
       xlab = "mean coverage", col = colb,
       breaks = seq(0, max(snp) + bz, bz))
  abline(v = thresholds, col = "red", lty = 2, lwd = 2)
  dev.off()
  }
  cat("Filter based on mean coverage (sum of AvgCountSnp and AvgCountRef)\n",
      "using thresholds of:", thresholds, "x", ", which correspond to",
      median_dev,
      "deviations from the median.",
      "\nInitial number of loci:", adegenet::nLoc(genlight),
      "\nFinal number of loci:", adegenet::nLoc(output_gen),
      "\nLoci removed:", adegenet::nLoc(genlight) - adegenet::nLoc(output_gen),
      "\n")
  return(output_gen)
  }
#takes genligth, plots individual missingness and returns same genlight
plot_missingness_individual <- function(genlight, plot_name ){
  #vector with individual missingness
  source("code/functions/ind_miss.r")
  vector_missingness <- ind_miss(genlight)
  #plot
  colb <- adjustcolor("orange", alpha.f = 0.5)
  bz <- 20
  pdf(file.path(getwd(), paste0(plot_name, ".pdf")), width = 5, height = 4)
  hist(vector_missingness, main = "",
       xlab = "Missing data per individual", col = colb,
       breaks = bz, xlim = c(0, 1))
  abline(v = ind_miss_thresh, col = "red", lty = 2, lwd = 2)
  dev.off()
  genlight
}
#takes genligth, plots locus missingness and returns same genlight
plot_missingness_locus <- function(genlight, plot_name ){
  #vector with individual missingness
  source("code/functions/ind_miss.r")
  vector_missingness <- 1 - adegenet::other(genlight)$loc.metrics$CallRate
  #plot
  colb <- adjustcolor("orange", alpha.f = 0.5)
  bz <- 20
  pdf(file.path(getwd(), paste0(plot_name, ".pdf")), width = 5, height = 4)
  hist(vector_missingness, main = "",
       xlab = "Missing data per locus", col = colb,
       breaks = bz, xlim = c(0, 1))
  abline(v = 1 - locus_callrate_threshold, col = "red", lty = 2, lwd = 2)
  dev.off()
  genlight
}
