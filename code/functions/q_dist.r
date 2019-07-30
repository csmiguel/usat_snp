#it returns similarity between ancestries for q_lists (see description)
q_dist <- function(arranged_dist = arranged_dist, species = NULL){
  #arranged_dist is q_list with 1st level being datasets and 2nd level being Ks
  #species to match in arranged_dist names c("hyla", "pelo")
  # if species is kept NULL, then arranged_dist is not filtered
  if (is.null(species)){
    h <- arranged_dist
  }else if (!is.null(species)){
    h <- arranged_dist[grep(species, names(arranged_dist))]
  }
  hh <-
  seq_along(h[[1]]) %>%
    sapply(function(y){
      lapply(h, function(x) x[[y]]) %>%
        starmie::averagePairWiseSimilarityH()
    })
  names(hh) <-  names(h[[1]])
  hh
  }
