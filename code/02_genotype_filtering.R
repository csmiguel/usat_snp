###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: filter DARrTseq genotypes from Hyla and Pelobates
#DESCRIPTION:
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
#  REQUIRED FILES:
#   Description: genlight objects with raw genotypes for Hyla and Pelobates
#   Inpath: data/intermediate/raw_genotypes.rds
#  OUTPUT:
#    Description: raw genotypes in genlight format
#    Outpath: data/intermediate/filtered_genotypes.rds
###.............................................................................
library(dartR)
library(dplyr)
library(magrittr)

source("code/functions/gl_stats.r")
source("code/functions/ind_miss.r")
source("code/functions/remove_replica.r")
source("code/functions/coverage_filt.r")
source("code/parameters/gl_cleaning.r")

input_path <- "data/intermediate/raw_genotypes.rds"
gen <- readRDS(file = input_path)

#summary statistics for raw data
raw_stats <- geno_stats(genlist = gen, data = "raw")

##################
#Filtering
sink(file = "data/intermediate/filtering.log")
  #remove replica from Pelobates:
      #for Pelobates Ana included a blind
      #replica named 108541 from sample 106854.
gen$dart_pelo <- remove_replica(gen = gen$dart_pelo, "108541", "106854")
  #filter SNPs

dart_filt <- grep("dart", names(gen)) %>% #only for dart genotypes
  {genfiltnames <<- names(gen)[.]; .} %>%
  seq_along() %>%
  lapply(function(i){
    sp <- paste(names(gen)[[i]])
    cat(paste0("#log for ", sp), "\n")
    filt_gen <-
    plot_missingness_individual(gen[[i]],
      paste0("data/final/Individual_missingness", sp)) %>%
    dartR::gl.filter.callrate(method = "ind",
    threshold = 1 - ind_miss_thresh, v = 5) %>%
    dartR::gl.filter.repavg(repavg_threshold, v = 5) %>%
    #allele balance: doi/full/10.1111/mec.14792 and 10.1111/1755-0998.12997
    all_balance_filter(tresholds = balance_tresholds, include_plot = T,
      plot_name = paste0("data/final/Allele_balance_filt_", sp)) %>%
    #filter on coverage: x deviations from the median
    total_coverage_filtering(median_dev = max_coverage,
      include_plot = T, lower_end = F,
      plot_name = paste0("data/final/Coverage_filtering", sp)) %>%
      plot_missingness_locus(
        paste0("data/final/Locus_missingness", sp)) %>%
    dartR::gl.filter.callrate(method = "loc",
      threshold = locus_callrate_threshold, v = 5, recalc = T) %>%
    dartR::gl.filter.secondaries(method = sec_method, v = 5) %>%
    dartR::gl.filter.maf(threshold = min_maf, v = 5) %>%
    {real_pop <<- pop(.); .} %>%
    #assume all the samples belong to the same large pop.
    `slot<-`("pop", value = factor(rep("all", length(.@ind.names)))) %>%
    dartR::gl.filter.hwe(alpha = hwe_alpha) %>%
    `slot<-`("pop", value = real_pop)
  }
)
sink()
names(dart_filt) <- genfiltnames

#statistics on filtered genotypes
filt_stats <- geno_stats(genlist = dart_filt, data = "filt")

#save
saveRDS(dart_filt, paste0("data/intermediate/dart_filt.rds"))
saveRDS(filt_stats, paste0("data/intermediate/filt_dartseq_stats.rds"))
saveRDS(raw_stats, paste0("data/intermediate/raw_dartseq_stats.rds"))
