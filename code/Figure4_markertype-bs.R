###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: Plot results from bootstrap analysis
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
#plot results from bootstrap analysis
bstree <- readRDS("data/intermediate/bs_treemetrics.rds")
# plot parameters
x_lab <- c("Distance to parent node", "Distance from root")
max_dnr <- sapply(bstree, function(x) x$dnr) %>% unlist %>% max
max_dnp <- sapply(bstree, function(x) x$dnp) %>% unlist %>% max
x_lim <- matrix(c(0, max_dnp, 0, max_dnr), ncol = 2)

#rename datasets
bstree
names(bstree) <-
  plyr::mapvalues(
    x = names(bstree),
    from = names(bstree),
    to = c("SNPs H. molleri", "SNPs P. cultripes",
    "microtellites H. molleri", "microtellites P. cultripes"))
#produce plots
for (i in 1:2){
  pdf(file = paste0("data/final/",
  gsub(" ", "_", x_lab[i]), ".pdf"),
   height = 7, width = 8.8)
  par(mfrow = c(2, 2), las = 1)
    for (p in 1:length(bstree)){
      hh <- bstree[[p]]
      plot(as.numeric(as.matrix(hh[, i])),
          hh$bs_tm, pch = 20,
          xlab = x_lab[i],
          ylab = "Bootstrap support",
          main = names(bstree)[p],
          xlim = x_lim[, i],
          ylim = c(0, 1))
      }
  dev.off()
}

#plot 2
plot_data <-
seq_along(bstree) %>%
  lapply(function(x){
  bstree[[x]] %>% mutate(dataset = as.factor(names(bstree)[x]))
    }) %>% {do.call(rbind, .)} %>%
  dplyr::mutate(species = as.character(dataset) %>%
                stringr::str_remove("[aA-zA]++ ")) %>%
  dplyr::mutate(marker = as.character(dataset) %>%
                stringr::str_remove(" .*$")) %>%
  {reshape2::melt(., id.vars = names(.)[-grep("dnr|dnp", names(.))])} %>%
    tibble::as_tibble() %>%
  dplyr::select(-dataset) %>%
  dplyr::rename(dist = variable)

library(cowplot)
library(ggplot2)

labels_distance <- c(dnp = "Dnp", dnr = "Dnr")
p <-
ggplot(data = plot_data) +
  geom_point(mapping = aes(x = value, y = bs_tm,
    colour = marker, shape = marker)) +
  facet_grid(species ~ dist, scales = "free_x",
             labeller = labeller(dist = labels_distance)) +
  labs(y = "Bootstrap support", x = "") +
  theme(
    strip.text.y = element_text(size = 10, face = "italic"),
    strip.text.x = element_text(size = 10)
  ) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  scale_shape_manual(values = c(5, 16)) +
  theme(legend.position = c(0.35, 0.3),
        legend.background = element_rect(size = 0.2, linetype = "solid",
                                 colour = "black"),
        legend.text = element_text(colour = "black", size = 10),
        legend.title =  element_blank())

p2 <- add_sub(p, label = c("Dnp"), x = 0.25, y = 3, vjust = 1, size = 10,
  vpadding = grid::unit(0, "lines"))
p3 <- add_sub(p2, label = c("Dnr"), x = 0.75, y = 3, vjust = 0, size = 10,
  vpadding = grid::unit(0, "lines"))
ggdraw(p3)
ggsave("data/final/Figure4_ggplot.pdf", height = 6, width = 8)
