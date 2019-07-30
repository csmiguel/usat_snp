###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: create podtopop populations_file for CLUMPAK with popID from Table1.rds
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(dplyr)
library(dartR)

input_path1 <- "data/intermediate/gen_consolidated_filtered.rds"
input_path2 <- "data/intermediate/table1.rds"
gen <- readRDS(file = input_path1)
meta <- readRDS(file = input_path2)

#structure runs datasets
str_files <- dir(path = "data/final", pattern = "run_")

#names from gen object that correspond to the dataset
gen_names <-
  str_files %>%
  replace(grep("pelo_[0-9]", str_files), "run_dart_pelo") %>%
  replace(grep("hyla_[0-9]", str_files), "run_dart_hyla") %>%
  gsub(pattern = "run_", replacement = "")

index_pop <- cbind(str_files, gen_names)

#for each dataset create a pop file
seq_along(str_files) %>%
lapply(function(x){
  #select element from list that matches names
  subset_gen <- gen[grep(index_pop[x, 2], names(gen))][[1]]
  #create table with IDs
  subset_gen$other$metadata %>%
    dplyr::left_join(meta, by = c("locality" = "locality")) %>%
    as.data.frame %>%
    dplyr::select(ID) %>%
    write.table(file = file.path("data/intermediate/clumpak",
                         paste0(index_pop[x, 1], "_populations_file")),
                col.names = F, quote = F, row.names = F)
})
