gl.stats <- function(genlight){
  x <- genlight
  firstc <- which(names(x@other$loc.metrics) == "CallRate")
  lastc <- ncol(x@other$loc.metrics)
  metrics <<- x@other$loc.metrics[, c(firstc:lastc)]
  mean1 <- apply(metrics, 2, mean)
  q1 <- apply(metrics, 2, quantile)
  sd1 <- apply(metrics, 2, sd)
  all1 <- cbind(mean1, t(q1), sd1)
  colnames(all1) <- c("mean", dimnames(q1)[[1]], "sd")
  return(all1)
}
