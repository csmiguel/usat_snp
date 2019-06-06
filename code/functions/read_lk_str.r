#read structure output for doing trace plots
run <- function(method = c("k1", "all")){
if (method == "k1")
r <- system("find data/final/runK1* -name *stlog",
  intern = T)
if (method == "all")
r <- system("find data/final/run_* -name *stlog",
  intern = T)
r
}

read_lk <- function(runs = runs){
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
  #change cols to numeric
  for (i in 2:ncol(p)){
    p[, i] <- as.numeric(as.character(p[, i]))
  }; rm(i)
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

plot_convergence <- function(method = c("with_burnin", "no_burnin"),
p = NULL, burnin = NULL){
  #p is an ouput from read_lk
  #burnin is the burnin from the parameter file
  it <- as.numeric(names(p)[3:length(names(p))])
  for (datset in 1:nlevels(p$dataset)){#for each dataset
    levelp <- levels(p$dataset)[datset] #name of dataset
    pd <- filter(p, dataset == levelp)#filter dataset
    pdf(paste0("data/final/convergence_", levelp, "_", method, ".pdf"))
    par(mfrow = n2mfrow(length(unique(pd$k))))
    for (K in min(pd$k):max(pd$k)){
      pdk <- filter(pd, k == K) %>% select(-dataset, -k)
      x_lim <- c(1, max(it))
      y_lim <- range(pdk)
      if(method == "no_burnin"){
        burnin_cols <- c(1: which(as.numeric(names(p)) == burnin))
        pdk <- filter(pd, k == K) %>% select(-burnin_cols)
        x_lim <- c(burnin, max(it))
        y_lim <- range(pdk)
      }
      c(1, max(it))
      plot(1, type = "n", xlim = x_lim, ylim = y_lim,
           main = paste0("K = ", K), ylab = "Ln Like",
           xlab = "Iterations")
      abline(v = burnin, col = "red", lty = 2)
      for (lin in 1:nrow(pdk)){
        lines(x = as.numeric(names(pdk)),  y = pdk[lin, ])
      }
    }
    dev.off()
  }
}


rhat <- function(p = NULL, burnin = NULL){
  it <- as.numeric(names(p)[3:length(names(p))])
  thinning_number <- it[2] - it[1] #frequency of saving steps
  #starting iteration after burnin
  start_itertion_number <- burnin + thinning_number
  burnin_cols <- c(1: which(as.numeric(names(p)) == burnin))
  p %>% select(-burnin_cols[-c(1, 2)]) %>%
  dplyr::mutate_if(is.factor, as.character) %>%
  plyr::dlply(c("dataset", "k"), function(x) {
  y <- t(as.matrix(x[, -c(1, 2), drop = FALSE]))
  y <- as.numeric(y)
  y <- as.list(as.data.frame(matrix(y, nrow = ncol(x) - 2)))
  l <- lapply(y, function(z) {
    coda::mcmc(z, start = start_itertion_number,
               thin = thinning_number)
  })
  gr <- coda::gelman.diag(coda::mcmc.list(l), autoburnin = FALSE)
  data.frame(id = as.character(x$dataset[1]),
             k = as.numeric(as.character(x$k[1])),
             rhat = gr$psrf[[1]], upper_ci =  gr$psrf[[2]])
}
) %>% plyr::rbind.fill()
}
