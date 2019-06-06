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
