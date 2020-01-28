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
               aes(x = longitude, y = latitude, color = sMLH, size = 1.8)) +
    scale_colour_gradientn(colours = c("blue", "white", "red"),
                         limits = c(0, 2)) +
    scale_fill_gradient(low = "gray90", high = "black") +
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

plot_res <-
  function(species = c("H. molleri", "P. cultripes"), plot_title = NULL){
    sp <-
      plyr::mapvalues(species, c("H. molleri", "P. cultripes"),
                  c("hyla", "pelo"))
  ggplot() +
    geom_tile(data = dplyr::filter(mapr_df, !is.na(value)),
              aes(x = x, y = y, fill = value)) +
    geom_path(data = dplyr::filter(distr, id == sp),
              aes(x = long, y = lat, group = group)) +
    geom_point(data = dplyr::filter(plot_residuals, dataset == species),
               aes(x = longitude, y = latitude, color = res, size = 2)) +
    scale_colour_gradientn(colours = c("darkblue", colm[2], "white",
                                      colm[1], "darkgreen"),
          limits = c(- max(plot_residuals$res), max(plot_residuals$res))) +
    scale_fill_gradient(low = "gray90", high = "black") +
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
