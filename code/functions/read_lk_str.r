#read structure output for doing trace plots
run <- function(method = c("k1", "all")){
if (method == "k1")
r <- system("find data/final/runK1* -name *stlog | grep -v usat",
  intern = T)
if (method == "all")
r <- system("find data/final/run_* -name *stlog | grep -v usat",
  intern = T)
r
}

read_lk <- function(){
  p <- runs %>%
    sapply(function(x){
    dataset <- gsub("^.*run[_]?(.*)/.*$", "\\1", x) %>% as.character()
    k <- gsub("^.*K([1-9]).*$", "\\1", x)
    lk <- system(paste0("cat ", x, " | grep [0,5][0]:"), intern = T) %>%
      #.[-grep("--", .)] %>%
      gsub("^.*(-\\d+).*$", "\\1", .) %>%
      c(dataset, k, .)
  }) %>% t %>% as.data.frame
}

get_lambda <- function(){
  runs %>%
    sapply(function(x){
  l <- system(paste0(
    "cat ", x, " | grep '^Mean value of lambda'| sed 's/Mean.*= //g'"),
    intern = T)
  dataset <- gsub("^.*run[_]?(.*)/.*$", "\\1", x) %>% as.character()
  c(dataset = dataset, lambda = l)
}) %>% t %>% dplyr::as_tibble() %>%
  dplyr::mutate(lambda = as.numeric(lambda)) %>%
  plyr::ddply("dataset", function(x) mean(x$lambda) %>% round(2))
  }
