source("code/parameters/structure.r")

#edit mainparams file for each STRUCTURE run

edit_mainparams <- function(gl = gl, name = name,
  str_file = str_file, dir_str = dir_str){
  param <- readLines("data/raw/mainparams")
  param[16] <- sub("xx",  10, param[16])#max pops (K)
  param[27] <- sub("xx",  adegenet::nInd(gl), param[27])#number of individual
  param[28] <- sub("xx",  adegenet::nLoc(gl), param[28])#number of loci
  param[22] <- sub("xx",  str_file, param[22])#infile.
  #structure_threader will take output and input infiles from -o and -i, but
  #not from the params files
    output <- file.path("data/final", paste0("str_", name, "_run"))
  param[23] <- sub("xx",  output, param[23])#outfile
  if (grepl("dart*", name) | is.numeric(name)){
    param[17] <- sub("xx", burnin_snp, param[17])#burnin length
    param[18] <- sub("xx", run_snp, param[18])#run length
  } else if (grepl("usat*", name)){
    param[17] <- sub("xx", burnin_usat, param[17])#burnin length
    param[18] <- sub("xx", run_usat, param[18])#run length
  }
  assertthat::assert_that(all(grepl("xx", param) == F),
    msg = "wrong replacement in mainparams file")
  mainparams_path <- file.path(dir_str, "mainparams")
  writeLines(param, mainparams_path)
}

edit_extraparams <- function(gl = gl, name = name,
  str_file = str_file, dir_str = dir_str){
  param <- readLines("data/raw/extraparams")
  if (grepl("dart*", name) | is.numeric(name)){
    param[93] <- sub("xx", run_snp / 200, param[93])#burnin length
  } else if (grepl("usat*", name)){
    param[93] <- sub("xx", run_usat / 200, param[93])#burnin length
  }
  assertthat::assert_that(all(grepl("xx", param) == F),
    msg = "wrong replacement in extraparams file")
  mainparams_path <- file.path(dir_str, "extraparams")
  writeLines(param, mainparams_path)
}
