###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2019
###.............................................................................
#GOAL: Plot convergence only for SNP subsets.
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................

rhat_results <- readRDS("data/intermediate/rhat.rds")
source("code/functions/read_lk_str.r")
source("code/parameters/structure.r")


library(ggplot2)
library(dplyr)

#Plot convergence for marker subsets of Pelobates
pelo_nloc <-
    readRDS("data/intermediate/gen_consolidated_filtered.rds")$dart_pelo %>%
     adegenet::nLoc()

pdf("data/final/rhat_subsampling.pdf", height = 4.5, width = 6)
plot(1, type = "n",
     xlim = range(rhat_results$dart$k),
     ylim = log10(range(rhat_results$dart$rhat)),
     ylab = "log10(rhat)",
     xlab = "K")
cols <- colorRampPalette(c("blue", "red"))
rhat_results$dart %>%
  filter(id != "dart_hyla") %>%
  mutate(id = as.character(id)) %>%
  mutate(id = gsub(pattern = "dart_pelo", pelo_nloc, id) %>% as.numeric()) %>%
  {ids <<- select(., id) %>% unique() %>% unlist %>% sort; .} %>%
  plyr::ddply("id", function(x) lines(x$k, log10(x$rhat),
                                      lwd = 2,
                                      type = "o",
                                      col = plyr::mapvalues(x$id,
                                            from = ids,
                                            to = cols(length(ids)))))
legend("topleft", legend = ids,
  col = cols(length(ids)), pch = 1, bty = "n", lwd = 2)
dev.off()
