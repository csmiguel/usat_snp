###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: export data into structure format
#DESCRIPTION: export data into structure format and create params files
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
#  REQUIRED FILES:
#   Description:
#   Inpath:
#  OUTPUT:
#    Description: raw genotypes plus metadata for dartseq samples
#    Outpath:
#  DEPENDENCIES:
###.............................................................................
library(dartR)
library(dplyr)

input_path <- paste0("data/intermediate/filt_genotypes.rds")
gen <- readRDS(file = input_path)
source("code/functions/edit_mainparams.r")

#warning dartR::gi2gl %>% gl2structure does not format format properly
#structure files. But, hierfstat::write.struct has a bug so that it does not
#print individual names in the ilab argument. I overcame this by modifiying the
#dataframe with genotypes:
names(gen) %>% seq_along() %>%
  sapply(function(x){
    dir_str <- file.path("data/intermediate", paste0("str_", names(gen)[x]))
    dir.create(dir_str, showWarnings = T)
    file.copy("data/raw/extraparams", file.path(dir_str, "extraparams"),
      overwrite = T) #same for all runs
    file.copy("data/raw/mainparams", file.path(dir_str, "mainparams"),
      overwrite = T) #edit for each run
    str_file <- file.path(dir_str, paste0(names(gen)[x], ".str"))
    edit_mainparams(gl = gen[[x]], name = names(gen)[x],
      str_file = str_file, dir_str = dir_str)
    if (class(gen[[x]]) == "genind"){
      hierfstat::genind2hierfstat(gen[[x]]) %>%
      dplyr::mutate(pop = as.factor(indNames(gen[[x]]))) %>%
      hierfstat::write.struct(fname = str_file)
      } else if (class(gen[[x]]) == "genlight") {
      dartR::gl2structure(gen[[x]], outfile = str_file)
    }
  })


#   cp extraparams to all subfolders
# read str stlog  cat K6_rep4.stlog | grep [00]: > lk.txt
# h <- read.table("data/intermediate/lk.txt", skip =1)
# plot(h$X0.068, type = "l")
