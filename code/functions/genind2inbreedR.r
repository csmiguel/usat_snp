#convert microsatellite genotypes in genind format to inbreedR
usatgenind2inbreedR <- function(x){
assertthat::assert_that(inherits(x, "genind"))
  l <-
    seq_along(x@loc.n.all) %>%
    sapply(function(y){
      lc <- x@loc.n.all
      x@tab[, (cumsum(lc)[y] - (lc[y] - 1)):cumsum(lc)[y]] %>%
        apply(1, function(zz) sort(zz, decreasing = T)[1])
    })
    l[l == 2] <- 0
    assertthat::assert_that(check_data2(l),
              msg = "incorrect inbreedR format")
    cat("genotypes formatted correctly for inbredR")
    return(l)
}

check_data2 <- function (genotypes, num_ind = NULL, num_loci = NULL){
  genotypes <- data.table::as.data.table(genotypes)
  vals <- unique(unlist(lapply(genotypes, unique)))
  if (length(unique(unlist(vals))) > 3) {
    stop("The data contains more than 3 elements (1, 0, missing)")
  }
  if (!all( (unique(unlist(vals))) %in% c(NA, 0, 1))) {
    stop("incorrect coding of matrix")
  }
  if (!is.null(num_ind)) {
    if (num_ind != nrow(genotypes)) {
      stop("Number of rows is unequal to the number of individuals")
    }
  }
  if (!is.null(num_loci)) {
    if (num_loci != ncol(genotypes)) {
      stop("Number of columns is unequal to the number of loci")
    }
  }
  return(TRUE)
}
