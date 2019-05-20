source("code/parameters/structure.r")
options(scipen = 10)
#edit mainparams file for each STRUCTURE run

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
