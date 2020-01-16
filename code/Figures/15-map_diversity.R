###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2019
###.............................................................................
#GOAL: Plot trees with heterozygosity
#PROJECT: usat_snp (https://github.com/csmiguel/usat_snp)
###.............................................................................
library(dartR)
library(dplyr)
library(maps)
library(ggplot2)
library(sf)
library(raster)

#load functions
source("code/functions/centroid_sf.r")
source("code/functions/map_diversity.r")
source("code/parameters/plotting_par.r")

#load data
#genotypes
gen <- readRDS("data/intermediate/gen_consolidated_filtered.rds")
#median sMLH
median_het <- readRDS("data/intermediate/median_het.rds")
#raster
p_raster <- "data/raw/raster_iberian_peninsula.grd"
mapr <- raster::raster(p_raster)
#IUCN distribution. Load shapefile and convert to dataframe
p_pelo <- "data/raw/redlist_species_data_44f009a9-b3b2-4608-a818-24c3dd714900"
p_hyla <- "data/raw/redlist_species_data_eff4cf11-63bb-4131-bba3-104c3358697b"
shp_pelo <-
  rgdal::readOGR(dsn = p_pelo, layer = "data_0") %>%
  sp::spTransform(mapr@crs) %>% raster::crop(extent(mapr)) %>%
  broom::tidy() %>% #convert to dataframe
  #I have to change the names of the factor levels so as
  # not to mix the polygons from both species.
  mutate(group = as.character(group) %>%
           stringr::str_replace("0.", "") %>% as.factor) %>%
  dplyr::mutate(id = "pelo")
shp_hyla <- rgdal::readOGR(dsn = p_hyla, layer = "data_0") %>%
  sp::spTransform(mapr@crs) %>% raster::crop(extent(mapr)) %>%
  broom::tidy() %>% #convert to dataframe
  dplyr::mutate(id = "hyla")

#bind both dataframes
distr <- rbind(shp_pelo, shp_hyla) %>%
  dplyr::mutate(id = as.factor(id))

#reformat raster to df
mapr_df <- cartomisc::gplot_data(mapr) #convert raster to dataframe

#format data from median heterozygosity
median_het_ggplot2 <-
  seq_along(median_het) %>%
    plyr::adply(1, function(x){
      datset_name <- names(median_het)[x]
      #heterozyosity
      het <-
        median_het[[datset_name]] %>%
        {data.frame(locality = names(.), sMLH = ., row.names = NULL)} %>%
        dplyr::mutate(dataset = names(median_het)[x])
      #compute centroid for each locality
      meta <-
        centroid_sf(gen[[datset_name]]@other$metadata) %>%
        sfc_as_cols
      #joint sMLH and coordinates by locality
      dplyr::left_join(het, meta, by = "locality")
      }) %>%
  rbind %>%
  dplyr::select(-1) %>%
  dplyr::as_tibble()

#individual plots
p_dart_hyla <- plot_frog(species = "hyla", genotypes = "dart_hyla",
                         plot_title = "SNPs H. molleri")
p_usat_hyla <- plot_frog(species = "hyla", genotypes = "usat_hyla",
                         plot_title = "microsatellites H. molleri")
p_dart_pelo <- plot_frog(species = "pelo", genotypes = "dart_pelo",
                         plot_title = "SNPs P. cultripes")
p_usat_pelo <- plot_frog(species = "pelo", genotypes = "usat_pelo",
                         plot_title = "microsatellites P. cultripes")

#get legend
legend1 <- cowplot::get_legend(p_dart_hyla)

#compose plot
all <-
  cowplot::plot_grid(
    p_dart_hyla + theme(legend.position = "none"),
    p_usat_hyla + theme(legend.position = "none"),
    p_dart_pelo + theme(legend.position = "none"),
    p_usat_pelo + theme(legend.position = "none"),
    ncol = 2, nrow = 2, align = "hv",
    labels = "AUTO")

cowplot::plot_grid(all, legend1, ncol = 2, rel_widths = c(1, 0.1))
ggsave("data/final/map_diversity.pdf", width = 9, height = 6)

#divesity map with residuals from 1:1 correlations
#residuals
residuals <- readRDS("data/intermediate/sMLH_residuals.rds")

plot_residuals <-
  #temporal object with unique coordinates per locality
  median_het_ggplot2[-which(duplicated(median_het_ggplot2$locality)), ] %>%
  #add coordinates to residuals
  {dplyr::left_join(residuals, ., by = "locality")} %>%
  dplyr::select(species, locality, res, longitude, latitude) %>%
  dplyr::rename(dataset = species)

p_res_hyla <-
  plot_res(species = "H. molleri", plot_title = "residuals H. molleri")
p_res_pelo <-
  plot_res(species = "P. cultripes", plot_title = "residuals P. cultripes")

#get legend
legend_res <- cowplot::get_legend(p_res_hyla)

#compose plot
all_res <-
  cowplot::plot_grid(
    p_res_hyla + theme(legend.position = "none"),
    p_res_pelo + theme(legend.position = "none"),
    ncol = 2, nrow = 1, align = "hv",
    labels = "AUTO")

cowplot::plot_grid(all_res, legend_res, ncol = 2, rel_widths = c(1, 0.2))
ggsave("data/final/map_diversity_residuals.pdf", width = 10, height = 5)

##report diversity
#H. molleri
#mean ± SD northern lineage
median_het_ggplot2 %>% dplyr::filter(dataset == "dart_hyla" & latitude > 41) %>% pull(sMLH) %>% {sapply(list(mean, sd), do.call, args = list(.))}
#mean ± SD southern lineage
median_het_ggplot2 %>% dplyr::filter(dataset == "dart_hyla" & latitude < 41) %>% pull(sMLH) %>% {sapply(list(mean, sd), do.call, args = list(.))}
#localities with lowest values
median_het_ggplot2 %>% dplyr::filter(dataset == "dart_hyla") %>% arrange(sMLH)
#localities with highest values
median_het_ggplot2 %>% dplyr::filter(dataset == "dart_hyla") %>% arrange(desc(sMLH))

#H. molleri
#mean ± SD northern lineage
median_het_ggplot2 %>% dplyr::filter(dataset == "dart_pelo") %>% arrange(sMLH)
