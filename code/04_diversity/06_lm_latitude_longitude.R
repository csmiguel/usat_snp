###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: explore latitude and longitude effects on genetic diversity
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(dplyr)
library(sf)

#function sfc_as_cols to convert sf to data frame
source("code/functions/centroid_sf.r")
#load datasets for heterozygosity, admixture and latitude
# and format them so to allow left join.
#1. sMLH
het <-
  readRDS("data/intermediate/sMLH_reformatted.rds") %>%
  dplyr::select(-sMLH_SNPs_median, -sMLH_usats_median, -locality) %>%
  reshape2::melt(id = c("species", "sample_id")) %>%
  dplyr::mutate(marker = as.character(variable) %>%
                  stringr::str_remove("sMLH_")) %>%
  dplyr::select(-variable) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::rename(smlh = value) %>%
  dplyr::as_tibble()

#2. latitude
# load genotypes with metadata
gen <-
  readRDS("data/intermediate/gen_consolidated_filtered.rds")
# create data frame with metadata with latitude
meta <-
  names(gen) %>%
  lapply(function(x) {
    h <-
      gen[grepl(x, names(gen))][[1]]
    h <- sfc_as_cols(h$other$metadata) %>%
      dplyr::mutate(species = ifelse(grepl("pelo", x),
                              "P. cultripes", "H. molleri"),
                    marker = ifelse(grepl("dart", x), "SNPs", "usats"))
    }) %>%
  do.call(what = rbind)

#3. coefficient of admixture
# load clumpak
qclumpak <-
  readRDS("data/intermediate/clumpak_major.rds") %>%
  {.[grep(pattern = c("usat_pelo|dart_pelo|dart_hyla|usat_hyla"),
                x = names(.))]}
# function for CA
source("code/functions/coeff_admixture.r")
# q-matrices do not have the names of the individuals. I will retrieve that from
assertthat::assert_that(all(sort(names(qclumpak)) == sort(names(gen))))

admx <-
  names(qclumpak) %>%
  lapply(function(x) {
    #integrate ca across K2-K8 (mean)
    h <- qclumpak[[x]] %>%
      lapply(function(y) apply(y, 1, ca)) %>%
      do.call(what = cbind) %>%
      apply(1, mean, na.rm = T)
    data.frame(ca = h,
               sample_id = adegenet::indNames(gen[[x]]),
               datset = x)
    }) %>%
    do.call(what = rbind) %>%
    dplyr::mutate(species = ifelse(grepl("pelo", datset),
                            "P. cultripes", "H. molleri"),
                  marker = ifelse(grepl("dart", datset), "SNPs", "usats")) %>%
    dplyr::select(-datset) %>%
    dplyr::as_tibble()

# join latitude, heterozygosity and coefficient of admixture
data <-
  dplyr::left_join(meta, het,
    by = c("sample_id" = "sample_id", "marker" = "marker")) %>%
  dplyr::left_join(admx,
    by = c("sample_id" = "sample_id", "marker" = "marker")) %>%
  dplyr::select(-species.x, -species.y) %>%
  dplyr::mutate(marker = as.factor(marker),
                species = as.factor(species))

#Models
#P. cultripes
d_pelo <- dplyr::filter(data, species == "P. cultripes")
# microsats
d_pu <- filter(d_pelo, marker == "usats")
m1_pelo_usat <- lm(smlh ~ longitude + latitude, data = d_pu)
# SNPs
d_ps <- filter(d_pelo, marker == "SNPs")
m1_pelo_dart <- lm(smlh ~ longitude + latitude, data = d_ps)

#H. molleri
d_hyla <- dplyr::filter(data, species == "H. molleri")
# microsats
d_hu <- filter(d_hyla, marker == "usats")
m1_hyla_usat <- lm(smlh ~ longitude + latitude, data = d_hu)
# SNPs
d_hs <- filter(d_hyla, marker == "SNPs")
m1_hyla_dart <- lm(smlh ~ longitude + latitude, data = d_hs)

sink("data/final/lm-models-diversty.txt")
cat(
  paste("Linear models to evaluate effect of Lat/Lon on genetic diversity\n",
"\n#H. molleri: microsatellites\n"))
summary(m1_hyla_usat)
cat("\n#H. molleri: SNPs\n")
summary(m1_hyla_dart)
cat("\n#P. cultripes: microsatellites\n")
summary(m1_pelo_usat)
cat("\n#P. cultripes: SNPs\n")
summary(m1_pelo_dart)
sink()

#plot
plot(smlh ~ ca, data = data, col = species, pch = as.character(marker))
