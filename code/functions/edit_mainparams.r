#edit mainparams file for each STRUCTURE run
edit_mainparams <- function(gl = gl, name = name,
  str_file = str_file, dir_str = dir_str){
  source("code/parameters/structure.r")
  h <- c("dart_hyla", "dart_pelo", "usat_hyla", "usat_pelo")
  assertthat::assert_that(name %in% h, msg = "check name in edit_mainparams")
  param <- readLines("data/raw/mainparams")
  param[16] <- sub("xx",  10, param[16])#number of individual
  param[27] <- sub("xx",  adegenet::nInd(gl), param[27])#number of individual
  param[28] <- sub("xx",  adegenet::nLoc(gl), param[28])#number of loci
  param[22] <- sub("xx",  str_file, param[22])#infile
    output <- file.path("data/final", paste0("str_", name, "_run"))
  param[23] <- sub("xx",  output, param[23])#outfile
  if (grepl("dart*", name)){
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
