rhat_results <- readRDS("data/intermediate/rhat.rds")
conv <- readRDS("data/intermediate/prob_convergence.rds")
source("code/functions/read_lk_str.r")
source("code/parameters/structure.r")


library(ggplot2)
library(dplyr)

#Plot convergence
plot_convergence(method = "with_burnin", p = conv$usat, burnin = burnin_usat)
plot_convergence(method = "with_burnin", p = conv$dart, burnin = burnin_snp)

#Plot rhat
p1 <-
  ggplot(data = rhat_results[[1]], aes(x = k, y = rhat, color = id)) +
  geom_line() +
  geom_hline(yintercept = 1.05, col = "red", lty = 2) +
  coord_cartesian(ylim = c(0.9, 1.2)) +
  ggtitle(names(rhat_results)[1])
p2 <-
  ggplot(data = rhat_results[[2]], aes(x = k, y = rhat, color = id)) +
  geom_line() +
  geom_hline(yintercept = 1.05, col = "red", lty = 2) +
  coord_cartesian(ylim = c(0.9, 1.2)) +
  ggtitle(names(rhat_results)[2])

pdf(paste0("data/final/rhat.pdf"), height = 6, width = 8)
  gridExtra::grid.arrange(p1, p2)
dev.off()

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
