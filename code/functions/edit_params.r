source("code/parameters/structure.r")
options(scipen = 10)
#edit mainparams file for each STRUCTURE run
#create extraparams file for structure input
edit_mainparams <- function(gl = gl, name = name,
  str_file = str_file, dir_str = dir_str){
  param <- readLines("data/raw/mainparams")
  hh <- grep("define MAXPOPS ", param)
  param[hh] <- sub("xx",  10, param[hh])#max pops (K)
  hh <- grep("#define NUMINDS ", param)
  param[hh] <- sub("xx",  adegenet::nInd(gl), param[hh])#number of individual
  hh <- grep("#define NUMLOCI ", param)
  param[hh] <- sub("xx",  adegenet::nLoc(gl), param[hh])#number of loci
  hh <- grep("#define INFILE ", param)
  param[hh] <- sub("xx",  str_file, param[hh])#infile.
  #structure_threader overrides these param on top
  hh <- grep("#define OUTFILE ", param)
    output <- file.path("data/final", paste0("str_", name, "_run"))
  param[hh] <- sub("xx",  output, param[hh])#outfile
  hhburn <- grep("#define BURNIN ", param)
  hhreps <- grep("#define NUMREPS ", param)
  if (grepl("dart*", name) | is.numeric(name)){
    param[hhburn] <- sub("xx", burnin_snp, param[hhburn])#burnin length
    param[hhreps] <- sub("xx", run_snp, param[hhreps])#run length
  } else if (grepl("usat*", name)){
    param[hhburn] <- sub("xx", burnin_usat, param[hhburn])#burnin length
    param[hhreps] <- sub("xx", run_usat, param[hhreps])#run length
  }
  assertthat::assert_that(all(grepl("xx", param) == F),
    msg = "wrong replacement in mainparams file")
  mainparams_path <- file.path(dir_str, "mainparams")
  writeLines(param, mainparams_path)
}

#create extraparams file for structure input
edit_extraparams <- function(gl = gl, name = name,
  str_file = str_file, dir_str = dir_str, mode = mode){
  param <- readLines("data/raw/extraparams")
  hh <- grep("#define UPDATEFREQ ", param)
  if (grepl("dart*", name) | is.numeric(name)){
    param[hh] <-
     sub("xx", run_snp / 200, param[hh])
  } else if (grepl("usat*", name)){
    param[hh] <- sub("xx", run_usat / 200, param[hh])#burnin length
  }
  hh <- grep("#define INFERLAMBDA ", param)
  if (mode == "K1"){
    param[hh] <- sub("xx", 1, param[hh])
  } else if (mode == "normal"){
    param[hh] <- sub("xx", 0, param[hh])
  }
  assertthat::assert_that(all(grepl("xx", param) == F),
    msg = "wrong replacement in extraparams file")
  mainparams_path <- file.path(dir_str, "extraparams")
  writeLines(param, mainparams_path)
}

#create structure inputs for dart and usat data
create_str_input <- function(mode = c("normal", "K1")){
  if (mode == "K1") gen <- gen[grepl("dart", names(gen))]
seq_along(gen) %>%
  sapply(function(x){
    if (mode == "K1") dirn <- "strK1_" else if (mode == "normal") dirn <- "str_"
    dir_str <- file.path("data/intermediate", paste0(dirn, names(gen)[x]))
    dir.create(dir_str, showWarnings = T)
    str_file <- file.path(dir_str, paste0(names(gen)[x], ".str"))
    edit_mainparams(gl = gen[[x]], name = names(gen)[x],
      str_file = str_file, dir_str = dir_str)
    edit_extraparams(gl = gen[[x]], name = names(gen)[x],
      str_file = str_file, dir_str = dir_str, mode = mode)
    if (class(gen[[x]]) == "genind"){
      hierfstat::genind2hierfstat(gen[[x]]) %>%
      dplyr::mutate(pop = as.factor(indNames(gen[[x]]))) %>%
      hierfstat::write.struct(fname = str_file)
      } else if (class(gen[[x]]) == "genlight") {
      dartR::gl2structure(gen[[x]], outfile = str_file)
    }
  })
}
#create structure inputs for downsampling of dart
create_str_input_subsets <- function(mode = c("normal", "K1")){
  pelo <- dartR::gi2gl(pelo)
1:length(s) %>%
sapply(function(x){
  int <- sample(1:adegenet::nLoc(pelo), s[x], replace = FALSE)
  s_pelo <- pelo[, int]
  if (mode == "K1") dirn <- "strK1_" else if (mode == "normal") dirn <- "str_"
  dir_str <- file.path("data/intermediate", paste0(dirn, s[x]))
  dir.create(dir_str, showWarnings = T)
  str_file <- file.path(dir_str, paste0(s[x], ".str"))
  edit_mainparams(gl = s_pelo, name = s[x],
      str_file = str_file, dir_str = dir_str)
  edit_extraparams(gl = s_pelo, name = s[x],
      str_file = str_file, dir_str = dir_str, mode = mode)
  if (class(s_pelo) == "genind"){
    hierfstat::genind2hierfstat(s_pelo) %>%
    dplyr::mutate(pop = as.factor(indNames(s_pelo))) %>%
    hierfstat::write.struct(fname = str_file)
    } else if (class(s_pelo) == "genlight") {
    dartR::gl2structure(s_pelo, outfile = str_file)
  }
  })
}

#notes:
#warning dartR::gi2gl %>% gl2structure does not format format properly
#structure files. But, hierfstat::write.struct has a bug so that it does not
#print individual names in the ilab argument. I overcame this by modifiying the
#dataframe with genotypes:
