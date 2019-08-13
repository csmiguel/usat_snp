###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: collect clumpak output for major cluster
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(dplyr)
library(adegenet)

#genotypes
gen <- readRDS("data/intermediate/shared_all_consolidated_filtered.rds")
#function to rearrange clumpak results
source("code/functions/rearrange_k.r")
#function to calculate centroids from localities
source("code/functions/centroid_sf.r")

#Collect Clumpak output:
datasets <-
  dir("data/intermediate/shared_clumpak", pattern = "run_", full.names = T) %>%
    {.[-grep(., pattern = ".zip|_file")]}

#1. create list with Major Cluster ancestries
qclumpak <-
  seq_along(datasets) %>%
  lapply(function(x){
    #kruns <- dir(path = datasets[x], pattern = "K=[0-9]")
    k <- list()
    for (i in 1:8){
      k[[i]] <-
        paste0(datasets[x], "/K=", i,
          "/MajorCluster/CLUMPP.files/ClumppIndFile.output") %>%
        readLines() %>%
        gsub(pattern = "^.*: ", replacement = "") %>%
        textConnection %>%
        read.table
      names(k)[i] <- paste0("K", i)
    }
    k
  })

names(qclumpak) <- datasets %>%
  gsub(pattern = "^.*\\/run_shared_", replacement = "")

#replace shared gen in all file
qq <- readRDS("data/intermediate/clumpak_major.rds")
seq_along(qclumpak) %>%
  lapply(function(x){
    qq[[match(names(qclumpak)[x], names(qq))]] <<- qclumpak[[x]]
    })
qclumpak <- qq
rm(qq)

#2. mean ancestries per population

#assert that metadata names are ok
assertthat::assert_that(seq_along(gen) %>% sapply(function(x){
  names(gen[[1]]@other$metadata) == c("sample_id", "locality",  "geometry")
}) %>%
all())
#vector to reorder datasets
dor <- c("dart_hyla", "usat_hyla", "dart_pelo", "usat_pelo")
gen <- gen[match(dor, names(gen))]

# remove subset of loci
qclumpak1 <- qclumpak %>%
  .[grep("^[a-z]", names(qclumpak))] %>%
  {.[match(dor, names(.))]}
qclumpak1 <- qclumpak1[!sapply(qclumpak1, is.null)]

# it returns a list with the same structure as qclumpak
#but with mean ancestries per location
mean_anc <-
  seq_along(qclumpak1) %>%
  lapply(function(x){
    qlist <- qclumpak1[[x]]
    namex <- names(qclumpak1)[x]
    meta <- gen[[namex]]@other$metadata
    h <-
      2:length(qlist) %>%
      lapply(function(y){
        split(qlist[[y]], meta$locality) %>%
          lapply(function(s) apply(s, 2, mean)) %>%
          do.call(what = rbind)
      })
    names(h) <- names(qclumpak1[[1]])[2:length(qlist)]
    h
  })
names(mean_anc) <- names(qclumpak1)

#3. rearrange ancestries between datasets (usat vs snps) for each species for
#   easiness of comparison when plotting.
hylak <- reorder_ancestries(species = "hyla", mean_anc = mean_anc)
pelok <- reorder_ancestries(species = "pelo", mean_anc = mean_anc)
mean_anc_arranged <- c(hylak, pelok)

#4. save objects
saveRDS(qclumpak, "data/intermediate/shared_clumpak_major.rds")
saveRDS(mean_anc, "data/intermediate/shared_population_ancestry.rds")
saveRDS(mean_anc_arranged,
  "data/intermediate/shared_population_ancestry_arranged.rds")
