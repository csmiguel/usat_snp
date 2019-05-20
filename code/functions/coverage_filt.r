#calculate difference in coverage. I will calculate the difference in coverage
# between alleles, units being number of times the most covered allele is
#different from the other. ie: 2.5 means that the allele with most converage
#has in average 2.5 times more coverage than the allele with less converage.
#This step was applied in 10.1111/1755-0998.12997
#Thank Elspeth in acknowledgements for sharing part of the code.
diff_coverage <- function(genlight_object){
  snp <- genlight_object$other$loc.metrics
  xx <- 1 / (1 - ( (abs(snp[, "AvgCountRef"] - snp[, "AvgCountSnp"])) /
            pmax(snp[, "AvgCountRef"], snp[, "AvgCountSnp"])))
  assertthat::assert_that(length(xx) == adegenet::nLoc(genlight_object))
          }
