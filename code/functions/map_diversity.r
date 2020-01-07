#plot function for map diversity

plot_frog <-
  function(species = c("hyla", "pelo"),
          genotypes = c("dart_hyla", "dart_pelo", "usat_hyla", "usat_pelo"),
          plot_title = NULL){
  ggplot() +
    geom_tile(data = dplyr::filter(mapr_df, !is.na(value)),
              aes(x = x, y = y, fill = value)) +
    geom_path(data = dplyr::filter(distr, id == species),
              aes(x = long, y = lat, group = group)) +
    geom_point(data = dplyr::filter(median_het_ggplot2, dataset == genotypes),
               aes(x = longitude, y = latitude, color = sMLH, size = 2)) +
    scale_colour_gradientn(colours = c("blue", "white", "red"),
                         limits = c(0, 2)) +
    scale_fill_gradient(low = "white", high = "black") +
    theme_classic() +
    guides(size = FALSE) +
    ylab(NULL) +
    xlab(NULL) +
    labs(fill = "Elevation (m)") +
    ggtitle(plot_title) +
    theme(
      plot.title = element_text(face = "plain", size = 12)
    ) +
    coord_quickmap()
}
