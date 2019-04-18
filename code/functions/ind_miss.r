#calculate missingness per individual from genlight object
ind_miss <- function(gl){
  apply(as.matrix(gl), 1, function(x){
  sum(is.na(x))
  }
  ) / adegenet::nLoc(gl)
}
