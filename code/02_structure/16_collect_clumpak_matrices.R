###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: collect clumpak matrices from Pairwise comparisons
# see Jakobsson and Rosenberg 2002 10.1093/bioinformatics/btm233
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

library(dplyr)

source("code/parameters/plotting_par.r")

krange <- c(2:8)
#1. create list Symmetric similarity coefficients SSC for each K.
# that is, for each K I will read the distance matrices returned by CLUMPAK.
qclumpak <-
  lapply(krange, function(k){ #for each K
    #read CLUMPP distance matrix
    distmatrix <-
    paste0("data/intermediate/clumpak_comparison/K=", k,
      "/CLUMPP.files/ClumppPairMatrix") %>%
    readLines() %>%
    textConnection %>%
    read.table
    #read list of structure runs used to compute the matrices.
    runnames <- paste0("data/intermediate/clumpak_comparison/K=", k,
      "/CLUMPP.files/FilesToIndex") %>%
      readLines() %>%
      {sapply(strsplit(., "\t"), `[`, 2)}
    #add names to matrix
    rownames(distmatrix) <- runnames
    colnames(distmatrix) <- runnames
    #convert distance matrix to 3-column matrix
    h <- reshape2::melt(as.matrix(distmatrix), varnames = c("row", "col"))
    hh <- h[as.numeric(h$row) > as.numeric(h$col), ]
    hh
    })
names(qclumpak) <- paste0("K", krange)
#create data for plotting
plot_data <-
  plyr::ldply(seq_along(qclumpak), function(i) {
    x <- qclumpak[[i]]
    x$cluster <- names(qclumpak)[i]
    x
  }) %>%
  as_tibble() %>%
#add column according to if comparison was usat-usat, snp-snp or snp-usat
    dplyr::mutate(comp = ifelse(grepl("usat", row) & grepl("usat", col), "usat",
                         ifelse(grepl("dart", row) & grepl("dart", col), "dart",
                          "both")))

#save data with Symmetric similarity coefficients ready to plot
saveRDS(plot_data, "data/intermediate/SSC.rds")
