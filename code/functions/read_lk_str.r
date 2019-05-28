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
      {iteration <<- gsub(":.*$", "\\1", .) %>% gsub(" ", "", .); .} %>%
      gsub("^.*(-\\d+).*$", "\\1", .) %>%
      c(dataset, k, .)
  }) %>% t %>% as.data.frame
  names(p) <- c("dataset", "k", iteration)
  p
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


plot_convervenge <- function(method = c("with_burnin", "no_burnin")){
  for (datset in 1:nlevels(p$dataset)){#for each dataset
    levelp <- levels(p$dataset)[datset] #name of dataset
    pd <- filter(p, dataset == levelp)#filter dataset
    pdf(paste0("data/final/convergence", levelp, "_", method, ".pdf"))
    par(mfrow = n2mfrow(length(unique(pd$k))))
    for (K in min(pd$k):max(pd$k)){
      pdk <- filter(pd, k == K) %>% select(-dataset, -k)
      x_lim <- c(1, max(it))
      y_lim <- range(pdk)
      if(method == "no_burnin"){
        pdk <- filter(pd, k == K) %>% select(-burnin_cols)
        x_lim <- c(burnin_snp, max(it))
        y_lim <- range(pdk)
      }
      c(1, max(it))
      plot(1, type = "n", xlim = x_lim, ylim = y_lim,
           main = paste0("K = ", K), ylab = "Ln Like",
           xlab = "Iterations")
      abline(v = burnin_snp, col = "red", lty = 2)
      for (lin in 1:nrow(pdk)){
        lines(x = as.numeric(names(pdk)),  y = pdk[lin, ])
      }
    }
    dev.off()
  }
}
