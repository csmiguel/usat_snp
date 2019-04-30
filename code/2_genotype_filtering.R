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

source("code/functions/gl_stats.r")
source("code/functions/ind_miss.r")
source("code/functions/remove_replica.r")

input_path <- paste0("data/intermediate/raw_genotypes.rds")
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
source("code/parameters/gl_cleaning.r")
gen_filt <- grep("dart", names(gen)) %>% #only for dart genotypes
  {genfiltnames <<- names(gen)[.]; .} %>%
  seq_along() %>%
  lapply(function(i){
    sp <- paste(names(gen)[[i]])
    cat(paste0("#log for ", sp), "\n")
    dartR::gl.filter.callrate(gen[[i]], method = "ind",
    threshold = 1 - ind_miss_thresh, v = 5) %>%
    dartR::gl.filter.repavg(repavg_threshold, v = 5) %>%
    dartR::gl.filter.callrate(method = "loc",
      threshold = locus_callrate_threshold, v = 5, recalc = T) %>%
    dartR::gl.filter.secondaries(method = sec_method, v = 5) %>%
    dartR::gl.filter.maf(threshold = min_maf, v = 5) %>%
    {real_pop <<- pop(.); .} %>%
    `slot<-`("pop", value = factor(rep("all", length(.@ind.names)))) %>%
    dartR::gl.filter.hwe(alpha = hwe_alpha) %>%
    `slot<-`("pop", value = real_pop)
  }
)
sink()
names(gen_filt) <- genfiltnames

#statistics on filtered genotypes
filt_stats <- geno_stats(genlist = gen_filt, data = "filt")

######################
#Raw data
saveRDS(raw_stats, paste0("data/intermediate/raw_dartseq_stats.rds"))
#Filtered data
  #create objects:
  for (i in which(!(names(gen) %in% names(gen_filt)))){
    gen_filt[[i]] <- gen[[i]]
    names(gen_filt)[i] <- names(gen)[i]
  } #add usat genotypes to gen_filt object
saveRDS(gen_filt, paste0("data/intermediate/filt_genotypes.rds"))
saveRDS(filt_stats, paste0("data/intermediate/filt_dartseq_stats.rds"))
