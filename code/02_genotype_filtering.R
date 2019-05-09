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
#remove individuals from metadata that were filtered in genotypes
  meta_path <- paste0("data/intermediate/metadata_dartseq.rds")
  metadata_dartseq <- readRDS(file = meta_path)
#assert_that all ids in genotypes are present in metadata
    assertthat::assert_that(all(indNames(gen_filt$dart_hyla) %in%
      metadata_dartseq$metadata_darthyla$sample_id))
    assertthat::assert_that(all(indNames(gen_filt$dart_pelo) %in%
      metadata_dartseq$metadata_dartpelo$sample_id))
    assertthat::assert_that(
    #assert_that all ids unique
    list(indNames(gen_filt$dart_hyla), indNames(gen_filt$dart_pelo),
            metadata_dartseq$metadata_dartpelo$sample_id,
            metadata_dartseq$metadata_darthyla$sample_id) %>%
      sapply(function(x) all(duplicated(x)) == F) %>% all()
    )
    #remove individuals from metadata
metadata_dartseq$metadata_darthyla %<>%
  filter(sample_id %in% indNames(gen_filt$dart_hyla))
metadata_dartseq$metadata_dartpelo %<>%
  filter(sample_id %in% indNames(gen_filt$dart_pelo))
    #create list with metadata for all samples
  metadata_dartseq[[3]] <- gen$usat_hyla@other
  metadata_dartseq[[4]] <- gen$usat_pelo@other
  names(metadata_dartseq)[c(3, 4)] <- c("usat_hyla", "usat_pelo")
assertthat::assert_that(length(metadata_dartseq) == length(gen),
  msg = "metadata_all has a different length than gen_filt")
#####################
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
saveRDS(metadata_dartseq, "data/intermediate/metadata_all.rds")
