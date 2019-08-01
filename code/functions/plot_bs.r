#plot bootstrap support against Dnr and Dns
plot_bs <- function(mode = c("color", "bw")){
  if (mode == "color") col1 <- rev(colm) else if (mode == "bw") col1 <- rev(colmbw)
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
    scale_color_manual(values = col1) +
    scale_shape_manual(values = c(16, 16)) +
    theme(legend.position = c(0.35, 0.3),
           legend.background = element_rect(size = 0.2, linetype = "solid",
                                            colour = "black"),
           legend.text = element_text(colour = "black", size = 10),
           legend.title =  element_blank())

  p <- cowplot::add_sub(p, label = c("Dnp"), x = 0.25, y = 3, vjust = 1,
    size = 12, vpadding = grid::unit(0, "lines")) %>%
    cowplot::add_sub(label = c("Dnr"), x = 0.75, y = 3, vjust = 0, size = 12,
                vpadding = grid::unit(0, "lines"))
  cowplot::ggdraw(p)
  ggsave(paste0("data/final/DnrDnp_ggplot_", mode, ".pdf"),
    height = 6, width = 8)
}
