loci.stats <- function(x){
  firstc <- which(names(x@other$loc.metrics) == "CallRate")
  lastc <- which(names(x@other$loc.metrics) == "RepAvg")
  metrics <- x@other$loc.metrics[, c(firstc:lastc)]
  mean1 <- apply(metrics, 2, mean)
  q1 <- apply(metrics, 2, quantile)
  sd1 <- apply(metrics, 2, sd)
  all1 <- cbind(mean1, t(q1), sd1)
  colnames(all1) <- c("mean", dimnames(q1)[[1]], "sd")
  return(all1)
}

geno_stats <- function(pattern = "dart", genlist, data = c("raw", "filt")){
  grep(pattern, names(genlist)) %>% #only for dart genotypes
  seq_along() %>%
  lapply(function(i){
    sp <- paste(names(genlist)[[i]])
    summary_metrics <- loci.stats(genlist[[i]]) #metrics per loci
    ind_miss_raw <- ind_miss(genlist[[i]]) #missingness per individual
      #plot individual missingness
    pdf(file = paste0("data/intermediate/", data, "_", sp,
    "_ind_missingness.pdf"), width = 6, height = 4)
          plot(sort(ind_miss_raw),
            main = paste0("Sorted per indiviudal missingness in ", data,
            " data from ", sp),
            ylab = "Proportion of missing loci")
    dev.off()
    list(genotypes = paste0(data, "_", sp), summary_metrics = summary_metrics,
      ind_miss_raw = ind_miss_raw)
    }
  )
}
