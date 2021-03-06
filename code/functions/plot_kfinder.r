library(dplyr)

#function to plot result from Kfinder saved in kfinder_ls.rds

plot_k <- function(klist = NULL, index = NULL, runs_pattern = NULL, mode = 1){
  #runs_pattern, is the pattern to grep from names
  #index is Pritchard, Evanno or Parsimony
  #klist is list with data frames of k statistics
  to_plot <- klist[grep(runs_pattern, names(klist))]
  if (mode == "subsampling"){
    to_plot <-
      to_plot[names(to_plot) %>%
                stringr::str_extract("[0-9]+") %>% as.numeric() %>% order]
  }
  y_lim <- to_plot %>% sapply(function(x) x[[index]]) %>%
    unlist %>% {suppressWarnings(as.numeric(.))} %>% range(na.rm = T)

  ylabs <- plyr::mapvalues(
    x = index, warn_missing = F,
    from = names(klist[[1]])[-1],
    to = c("Pr[X|K]", "Evanno ΔK", "Parsimony index")
  )

  plot(1, type = "n",
       xlim = range(to_plot[[1]]$K),
       ylim = y_lim,
       ylab = ylabs,
       xlab = "K")
  #cols <- colorRampPalette(c("#d73027", "#1a9850"))
  #order of colors matches order of list to plot
  assertthat::assert_that(grep(markers[1], names(to_plot)) == 1)
  for (i in seq_along(to_plot)){
    lines(to_plot[[i]]$K,
          suppressWarnings(as.numeric(as.character(to_plot[[i]][[index]]))),
          lwd = 2,
          lty = ltym[i],
          type = "o",
          col = colm[i],
          pch = pchm[i]
          )
  }
  legend("bottomright", legend = names(to_plot),
         col = colm, pch = pchm, lwd = 2, lty = ltym,
         bg = adjustcolor("white", alpha.f = 0.5))
}
