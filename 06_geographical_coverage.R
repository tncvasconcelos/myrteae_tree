# rm(list=ls())
library(ape)
library(phytools)
library(data.table)
library(maptools)
library(raster)

# setwd("~/Desktop/myrteae_traits")

#########################
organize.bubble.plot <- function(trait_table, reference_table, all_vars, twgd_data) {
  tmp_reference_table <- subset(reference_table, reference_table$wcvp_name %in% unique(trait_table$species))
  wcvp_subset <- subset(all_vars, all_vars$taxon_name %in% tmp_reference_table$wcvp_name)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$introduced==0)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$extinct==0)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$location_doubtful==0)
  
  focal_areas <- unique(wcvp_subset$area_code_l3)
  results <- matrix(nrow=0, ncol=5)
  for(i in 1:length(focal_areas)) {
    one_area <- focal_areas[i]
    one_subset <- subset(wcvp_subset, wcvp_subset$area_code_l3==one_area)
    sp_rich <- length(unique(one_subset$taxon_name))
    family_rich <- length(unique(one_subset$family))
    area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% one_area),]
    if(nrow(area_plus_buffer)>0) {
      centroids <- rgeos::gCentroid(area_plus_buffer, byid=TRUE)
      lon <- extent(centroids)[1]
      lat <- extent(centroids)[3]
      results <- rbind(results, cbind(sp_rich, family_rich, one_area, lon, lat))
    }
    cat(i, "\r")
  }
  results <- as.data.frame(results)
  results$sp_rich <- as.numeric(results$sp_rich)
  results$family_rich <- as.numeric(results$family_rich)
  results$lon <- as.numeric(results$lon)
  results$lat <- as.numeric(results$lat)
  return(results)
}


#########################
#########################
dist_sample <- read.table("WCVP/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("WCVP/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

#########################
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")
reference_table <- readRDS("myrteae_species_WCVP_GBIF.Rdata")
myrteae_genera <- read.csv("neotropical_myrteae_genera_WCVP.csv")
tree <- read.tree("mmc_target_common_June27_pruned_no_out.tre")

all_vars <- subset(all_vars, all_vars$family %in% "Myrtaceae")
all_vars <- subset(all_vars, all_vars$taxon_status %in% "Accepted")
all_vars <- subset(all_vars, all_vars$species != "")
all_vars <- subset(all_vars, all_vars$genus %in% myrteae_genera$neotropical_myrteae_genera)
all_species <- unique(all_vars$taxon_name)

#-----------------------------
path="WCVP/wgsrpd-master/level3/level3.shp"
twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))
#-----------------------------

all_areas <- unique(all_vars$area)
results_geo <- data.frame(area=all_areas, prop_sampled=NA)
for(i in 1:length(all_areas)) {
  one_subset <- all_vars[which(all_vars$area==all_areas[i]),]
  sampled <- gsub("\\."," ", tree$tip.label)
  prop_sampled <- length(which(one_subset$taxon_name %in% sampled)) / nrow(one_subset)
  results_geo[i,2] <- prop_sampled
}

result_sampled <- data.frame(species=all_species, prop_sampled=NA)
for(u in 1:length(all_species)) {
  sampled <- gsub("\\."," ", tree$tip.label)
  if(all_species[u] %in% sampled) {
    result_sampled[u,2] <- "sampled"
  } else {
    result_sampled[u,2] <- "not_sampled"
  }
}

reference_table <- readRDS("myrteae_species_WCVP_GBIF.Rdata")

organized_table_for_plot_total <- organize.bubble.plot(result_sampled, reference_table, all_vars, twgd_data)
organized_table_for_plot_sampled <- organize.bubble.plot(subset(result_sampled, result_sampled$prop_sampled=="sampled"), reference_table, all_vars, twgd_data)

nas <- rep(NA, nrow(organized_table_for_plot_sampled))
proportion_table <- data.frame(sp_rich_prop=nas, sp_rich_total=nas, one_area=nas, lon=nas, lat=nas)
for(i in 1:nrow(organized_table_for_plot_sampled)) {
  annual_sp_rich <- organized_table_for_plot_sampled$sp_rich[i]
  one_area <- organized_table_for_plot_sampled$one_area[i]
  total_sp_rich <- organized_table_for_plot_total$sp_rich[organized_table_for_plot_total$one_area == one_area]
  one_proportion <- round(annual_sp_rich / total_sp_rich, 3)
  proportion_table$sp_rich_prop[i] <- one_proportion
  proportion_table$sp_rich_total[i] <- total_sp_rich
  proportion_table$one_area[i] <- one_area
  proportion_table$lon[i] <- organized_table_for_plot_sampled$lon[i]
  proportion_table$lat[i] <- organized_table_for_plot_sampled$lat[i]
}

library(ggplot2)
library(maps)
library(ggthemes)
library(viridis)

twgd_data01 <- sf::st_as_sf(twgd_data)
twgd_data01 <- merge(twgd_data01, proportion_table, by.x="LEVEL3_COD", by.y="one_area")
twgd_data_americas <- subset(twgd_data01, twgd_data01$LEVEL1_COD%in%c(7,8))

tmp_map1 <- ggplot(data = twgd_data_americas) +
  geom_sf(aes(fill = sp_rich_prop)) +
  scale_fill_viridis_c(option = "viridis", alpha=0.8, direction=-1) +
  theme_classic() 

tmp_map2 <- ggplot(data = twgd_data_americas) +
  geom_sf(aes(fill = sp_rich_total)) +
  scale_fill_viridis_c(option = "viridis", alpha=0.8, direction=-1) +
  theme_classic() 

pdf("prop_species_sampled.pdf")
tmp_map1
tmp_map2
dev.off()
