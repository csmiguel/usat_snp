#plot bootstrap support against Dnr and Dns
plot_bs <- function(mode = c("color", "bw")){
  if (mode == "color"){
    col1 <- rev(colm)
    col2 <- col1
  } else if (mode == "bw"){
    col1 <- rev(colmbw)
    col2 <- c("darkgrey", "black")
  }
  p <-
    ggplot(plot_data) +
    geom_point(aes(x = value, y = bs_tm, colour = marker, shape = marker)) +
    geom_smooth(aes(x = value, y = bs_tm, colour = marker),
                method = "lm") +
    facet_grid(species ~ dist, scales = "free_x",
               labeller = labeller(dist = labels_distance,
                                   species = labels_species)) +
    labs(y = "Bootstrap support", x = "") +
    cowplot::theme_cowplot() +
    theme(strip.text.y = element_text(size = 10, face = "italic"),
          strip.text.x = element_text(size = 10)) +
    scale_color_manual(values = col1,
                       labels = c("SNPs", "microsatellites")) +
    scale_fill_manual(values = col2, guide = FALSE) +
    scale_shape_manual(values = c(16, 16), guide = FALSE) +
    theme(legend.position = c(.3, .6),
          legend.background = element_rect(size = 0.2,
                                           linetype = "solid",
                                           colour = "black"),
          legend.text = element_text(colour = "black", size = 9),
          legend.title =  element_blank(),
          strip.text.x = element_blank())

  p <- cowplot::add_sub(p, label = c("Dnp"), x = 0.25, y = 3, vjust = 1,
    size = 12, vpadding = grid::unit(0, "lines")) %>%
    cowplot::add_sub(label = c("Dnr"), x = 0.75, y = 3, vjust = 0, size = 12,
                vpadding = grid::unit(0, "lines"))
  cowplot::ggdraw(p)
  ggsave(paste0("data/final/DnrDnp_ggplot_", mode, ".pdf"),
    height = 6, width = 8)
}
