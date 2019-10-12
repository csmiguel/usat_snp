#coefficient of admixture
ca <- function(x){
  camin <- length(x) * (1 / length(x)) ^ 2
  #x is a numeric vector with invidual memberships across k clusters
  (sum(x ^ 2) - camin) / (1 - camin)
}
