#for 2 datasets it uses clumpp to align q-matrices for a list of K.
reorder_ancestries <- function(species = c("hyla", "pelo"),
            mean_anc = mean_anc, kmax = NULL){
# species, greps species name from object with mean ancestries per population
# mean_anc, is a list of lists: 1st level is dataset, 2nd level are K clusters
# kmax, max K to rearrange
  #structure of mean ancestries is a list of lists
  assertthat::assert_that(all(sapply(mean_anc, is.list)))
  #create object from mean acestries filtering by species
  anc <- mean_anc[grep(species, names(mean_anc))]
  #assert that only 2 datasets are compared
  assertthat::assert_that(length(anc) == 2)
  #assert that the snp data is in position one (does not get columsn changed)
  #and usat is in position 2 (its columns are rearranged by starmie::clumpp
  #according to snp data)
  assertthat::assert_that(grep("dart",
        names(anc)) == 1 & grep("usat", names(anc)) == 2)
  if (!is.null(kmax)){
    kvalues <- 1:(kmax - 1)
  }else if (is.null(kmax)){
    kvalues <- seq_along(anc[[1]])
  }

  #use clumpp from starmie to align clusters between usats and snps for each K:
  results_clumpp <-
    seq_along(kvalues) %>% #for each K
      lapply(function(x){
        #create a list of matrices to compare with same k for both datasets
        h <- list(anc[[1]][[x]], anc[[2]][[x]])
        names(h) <- names(anc)
        starmie::clumpp(h)
      })
  names(results_clumpp) <- names(anc[[1]])[kvalues] #remane with K clusters
  #rearrange results of clumpp into a list of lists: dataset::K
  #instead of K::dataset.
  rearranged_anc <- list()
  for (marker in seq_along(anc)){
    temp <- list()
    for (k in seq_along(results_clumpp)){
      temp[[k]] <- results_clumpp[[k]]$Q_list[[marker]]#rarrange matrices
    }
    names(temp) <- names(anc[[1]])[kvalues] #rename with K cluster names
    rearranged_anc[[marker]] <- temp
  }
  rm(temp, k, marker)
  names(rearranged_anc) <- names(anc) #rename with names of datasets

  rearranged_anc
}
