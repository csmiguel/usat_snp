###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: compute Simmetry Similarity Coefficient to compare ancestry matrices
# from STRUCTURE runs for Hyla molleri. The approach is fundamentally different
# as the one taken for P. cultripes. Since for Hyla we have different samples,
# we are obliged to compare localities instead of individuals. So first, I
# compute average ancestries per population for every STRUCTURE run. Second,
# I compute means per locality. Third, I CLUMPP results for every K using
# starmie built-in function. Fourth, I compute SSC  using starmie, which is
# supposed to compute the same metric as the one in CLUMPAK sofware.
# see Jakobsson and Rosenberg 2002 10.1093/bioinformatics/btm233
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(dplyr)
library(dartR)

#loead genotypes
gen <-
  readRDS("data/intermediate/gen_consolidated_filtered.rds") %>%
  {.[grep("hyla", names(.))]} %>%
  {setNames(., stringr::str_replace(names(.), "_hyla", ""))}

#read structure and compute ancestry means per population for every replicate run.
# the result is a list: marker(usat,dart):K(2:8):reps(1:10)
markers <- c("dart", "usat")
qlist <-
  lapply(markers, function(marker){
  #
  meta <- gen[[marker]]@other$metadata
  #list of files with structure output
  {ls_str <<- dir(paste0("data/final/run_", marker, "_hyla"),
    pattern = "*K[2-8]_rep.*_f", full.names = T); ls_str} %>%
  #vector with K's
  sapply(function(files) stringr::str_extract(files, "K[1-9]")) %>%
  {ks <<- unique(.); ks} %>%
  lapply(function(K){
    #for every K create a list with Qmatrices
    grep(pattern = K, ls_str, value = T) %>%
      lapply(function(reps){
        starmie::loadStructure(reps) %>% #load structure matrix
          starmie::getQ() %>%
          data.frame() %>%
          split(meta$locality) %>% #split samples per locality
          lapply(function(s) apply(s, 2, mean)) %>% #compute mean per locality
          do.call(what = rbind)
        }) %>% {setNames(., paste(marker, K, "rep", c(1:length(.)), sep = "_"))}
      }) %>% setNames(paste(marker, ks, sep = "_"))
  }) %>% setNames(markers)

#clumpp + cumpute Pairwise distances using SCC for Hyla
ssc_hyla <-
  seq_along(ks) %>%
  lapply(function(x){
    h <-
      #create list with combined usat snp K matrices for each K
      c(qlist[["dart"]][[x]], qlist[["usat"]][[x]]) %>%
      starmie::clumpp() %>% #CLUMPP all matrices
      {.$Q_list}
    #create dataframe with all the combinations of matrices to compare
    dfSCC <- names(h) %>% combn(2) %>% t %>%
      data.frame(stringsAsFactors = F)
    #fill in a 3rd column with the SCC
    dfSCC$SCC <-
      apply(dfSCC, 1, function(s){
        starmie::averagePairWiseSimilarityH(list(h[[s[1]]],h[[s[2]]]))
      })
    #add column according to if comparison was usat-usat, snp-snp or snp-usat
    dplyr::as_tibble(dfSCC) %>%
      dplyr::mutate(comp = ifelse(grepl("usat", X1) & grepl("usat", X2), "usat",
                           ifelse(grepl("dart", X1) & grepl("dart", X2), "dart",
                                         "both")))
  }) %>% setNames(ks)

#save results
saveRDS(ssc_hyla, "data/intermediate/ssc_hyla.rds")
