####.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: Plot map wih samples
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

library(dartR)
library(dplyr)
#which samples are shared between datasets:
gen <- readRDS("data/intermediate/gen_consolidated_filtered.rds")
samples <- readRDS("data/intermediate/table1.rds")

#Individuals per dataset
indHm <- adegenet::nInd(gen$usat_hyla)
indHs <- adegenet::nInd(gen$dart_hyla)
indPm <- adegenet::nInd(gen$usat_pelo)
indPs <- adegenet::nInd(gen$dart_pelo)
#Individuals per species
#H. molleri
indH <-
  c(adegenet::indNames(gen$dart_hyla), adegenet::indNames(gen$usat_hyla)) %>%
  unique() %>% length()

#P. cultripes
indP <-
  c(adegenet::indNames(gen$dart_pelo), adegenet::indNames(gen$usat_pelo)) %>%
  unique() %>% length()

#total individuals
totalind <- indP + indH

#localities for H. molleri
popH <- samples$samples %>% grep(pattern = "H.m.") %>% length
#localities for P. cultripes
popP <- samples$samples %>% grep(pattern = "P.c.") %>% length
#missing data microtellite
# H. molleri
missHm <- hierfstat::genind2hierfstat(gen$usat_hyla) %>%
  dplyr::select(-1) %>% {xx <<- length(as.matrix(.)); . } %>%
  is.na() %>% {sum(.) / xx * 100} %>% round(1)
# P. cultripes
missPm <- hierfstat::genind2hierfstat(gen$usat_pelo) %>%
  dplyr::select(-1) %>% {xx <<- length(as.matrix(.)); . } %>%
  is.na() %>% {sum(.) / xx * 100} %>% round(1)
#missing data SNPs
missHs <- as.matrix(gen$dart_hyla) %>%
  {xx <<- length(.); . } %>%
  is.na() %>% {sum(.) / xx * 100} %>% round(1)
missPs <- as.matrix(gen$dart_pelo) %>%
  {xx <<- length(.); . } %>%
  is.na() %>% {sum(.) / xx * 100} %>% round(1)
#Number of loci
# H. molleri
#  microtellites
locHm <- adegenet::nLoc(gen$usat_hyla)
#  SNPs
locHs <- adegenet::nLoc(gen$dart_hyla)
# P. cultripes
#  microtellites
locPm <- adegenet::nLoc(gen$usat_pelo)
#  SNPs
locPs <- adegenet::nLoc(gen$dart_pelo)

#raw raw_genotypes
genraw <- readRDS("data/intermediate/raw_genotypes.rds")
rawHm <- adegenet::nInd(genraw$usat_hyla)
rawHs <- adegenet::nInd(genraw$dart_hyla)
rawPm <- adegenet::nInd(genraw$usat_pelo)
rawPs <- adegenet::nInd(genraw$dart_pelo)


#shared samples between datasets

indHoverlap_p <-
  adegenet::indNames(gen$dart_hyla) %in% adegenet::indNames(gen$usat_hyla) %>%
  {indHoverlap_n <<- sum(.); .} %>%
  {indHoverlap_n / length(.)} %>%
  {round(. *100, 1)}

indPoverlap_p <-
  adegenet::indNames(gen$dart_pelo) %in% adegenet::indNames(gen$usat_pelo) %>%
  {indPoverlap_n <<- sum(.); .} %>%
  {indPoverlap_n / length(.)} %>%
  {round(. *100, 1)}

#print results
sink("data/final/results.txt")
cat(
  "Numbers of genotypes:\nAfter filtering, the dataset kept for analysis included a total", totalind, "individuals from", nrow(samples), "localities: ", indH," individuals of Hyla molleri (microtellites: n =" ,indHm,", ", locHm, "loci,", missHm," % missing data; SNPs: n = ", indHs, ",", locHs, " loci, ", missHs, " % missing data), from", popH, " localities, and ", indP," individuals of Pelobates cultripes (microtellites: n = ", indPm, ", ", locPm, "loci, ", missPm, "% missing data; SNPs: n = ", indPs, ",", locPs, "loci, ", missPs, " % missing data), from", popP, "localities\n")

cat("Shared samples between datasets:\nFor H. molleri,", indHoverlap_n," out of the", indHs, ", individuals (", indHoverlap_p ,"%) in the SNP dataset were also genotyped for microsatellites, whereas for P. cultripes, the,", indPoverlap_n, ", out of the", indPm, "individuals (", indPoverlap_p, "%) where also genotyped for microsatellites.")
sink()
