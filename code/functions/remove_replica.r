#function to remove replicated samples.
#it will remove the sample with most missing data from each replicate pair
#requires: genlight object (gen), a vector with original sample names (original),
# and a vector with names of replicated samples in the same order to the
#corresponding original names

remove_replica <- function(gen, original, replicated){
  miss <- ind_miss(gen)#calculte missingness per individual
  assertthat::assert_that(length(original) == length(replicated))
  assertthat::assert_that(all(original %in% indNames(gen)))
  assertthat::assert_that(all(replicated %in% indNames(gen)))
  for (i in 1:length(original)){
    max_miss <- sort(miss[which(original[i] == names(miss) |
    replicated[i] == names(miss))], decreasing = T)[1] %>% names()
    #get name of sample with highest missingness
    #drop that sample from genlight and recalculate metrics
    gen <- dartR::gl.drop.ind(gen, max_miss, recalc = F, mono.rm = TRUE, v = 5)
  }
  #gen <- dartR::gl.recalc.metrics(gen, v = 5)
  return(gen)
}
