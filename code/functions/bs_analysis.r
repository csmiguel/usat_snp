bs_analysis <- function(tr, bs, mode = c(1, 2)){
#calculate max length for each tree
  library(ggplot2)
  require(gridExtra)

  tr1 <- tr[[1]][[1]]
  tr2 <- tr[[1]][-1]
  bs <- bs / nboot

  maxd <- ape::cophenetic.phylo(tr1) %>% max()
  #standarize
  tr1$edge.length <- tr1$edge.length / maxd
  assertthat::assert_that(round(max(ape::cophenetic.phylo(tr1)), 3) == 1)

  #get minimium distance of each node to a tip
  dno <- ape::dist.nodes(tr1)[-c(1:Ntip(tr1)), 1:Ntip(tr1)] %>%
      apply(1, function(x) sort(x)[1]) #issue: check calculations
      #issue: it might be better to remove terminal edges
  #matrix with pairwise node distances excluding tips
  dnn <- ape::dist.nodes(tr1)[-c(1:Ntip(tr1)), -c(1:Ntip(tr1))] #issue: check calculations

  #check node names are sorted in increasing order
  assertthat::assert_that((as.numeric(rownames(dnn))[length(rownames(dnn))]  -
   as.numeric(rownames(dnn))[1] + 1) / length(rownames(dnn)) == 1)

  dnnmin <- numeric(nrow(dnn)) #issue: check calculations
  for (i in 1:ncol(dnn)){
    if (i == 1) dnnmin[i] <- dnn[2, 1:(2 - 1)] %>% min() else
    #get min egde distance from a given node to a more central one
    dnnmin[i] <- dnn[i, 1:(i - 1)] %>% min()
  }
  rm(i)
  #to remove 1st bs and correspoding values
  dnnmin <- dnnmin[-1]
  dno <- dno[-1]
  bs <- bs[-1]
  ####
  assertthat::assert_that(length(dnnmin) == length(dno))
  assertthat::assert_that(length(dno) == length(bs))
  #statistics
  dnnmin.glm <- summary(glm(bs~dnnmin, family = binomial()))
    dnnmin.glm.e <- dnnmin.glm$coefficients[2, 1] %>% round(3)
    dnnmin.glm.s <- dnnmin.glm$coefficients[2, 4] %>%
      formatC(format = "e", digits = 2)
  dno.glm <- summary(glm(bs~dno, family = binomial()))
    dno.glm.e <- dno.glm$coefficients[2, 1] %>% round(3)
    dno.glm.s <- dno.glm$coefficients[2, 4] %>%
      formatC(format = "e", digits = 2)
  #create color palette
  colorfun <- colorRamp(c("red", "blue"))
  pal <- colorfun(bs)
if (mode == 1){
  plot(dno, dnnmin, pch = 20,
      #col = grey(1 - (bs / max(bs))),
      col = rgb(pal[, 1], pal[, 2], pal[, 3], maxColorValue = 255),
      ylab = "Distance to the closest deeper node",
      xlab = "Depth in the tree",
      main = paste("BS across the tree ", names(tr)),
      xlim = c(0, .35), ylim = c(0, .13))
    text(0.1, 0.1, paste0("depth: p-value = ", dno.glm.s,
     "; estimate = ", dno.glm.e), cex = 0.7)
    text(0.1, 0.09, paste0("closest node: p-value = ", dnnmin.glm.s,
      "; estimate = ", dnnmin.glm.e), cex = 0.7)
    } else if (mode == 2){
      plot(dno, bs, pch = 20,
          xlab = "Depth in the tree",
          ylab = "Bootstrap support",
          main = names(tr),
          xlim = c(0, .35),
          ylim = c(0, 1))
        text(0.1, 0.1, paste0("depth: p-value = ", dno.glm.s,
         "; estimate = ", dno.glm.e), cex = 0.7)
    } else if (mode == 3){
      plot(dnnmin, bs, pch = 20,
          xlab = "Distance to the closest deeper node",
          ylab = "Bootstrap support",
          main = names(tr),
          xlim = c(0, .14),
          ylim = c(0, 1))
        text(0.1, 0.09, paste0("closest node: p-value = ", dnnmin.glm.s,
          "; estimate = ", dnnmin.glm.e), cex = 0.7)
    }
  }
