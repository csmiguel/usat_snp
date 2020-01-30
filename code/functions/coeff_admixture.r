#coefficient of admixture
ca <- function(x){
  #x is a numeric vector with invidual memberships across k clusters.
  # e.g.: x <- c(0.2, 0.3, 0.5)
  camin <- length(x) * (1 / length(x)) ^ 2
  (sum(x ^ 2) - camin) / (1 - camin)
}
