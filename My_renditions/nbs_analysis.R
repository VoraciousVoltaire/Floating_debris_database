
# Loading relevant packages
library(sf)
library(terra)


# Working with non-breeding season's data----

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/input_data/")
new_data_1 <- readRDS("test_2colonies.rds")
new_data_2 <- readRDS("test_2colonies_individ_info.rds")
indiv_merged_df <- merge(new_data_1, new_data_2, by = "individ_id")
relevant_new_data_1 <- dplyr::select(indiv_merged_df, individ_id, timestamp, lon, lat, loc_type, colony)
sf_relevant <- st_as_sf(relevant_new_data_1[,c(3,4,6)], coords = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Setting a new crs based on median positions of geolocations----

mean_loc <- geosphere::geomean(cbind(relevant_new_data_1$lon, relevant_new_data_1$lat))
proj.laea = paste("+proj=laea +lat_0=",mean_loc[2], " +lon_0=",mean_loc[2]," +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km ", sep="")
sf_relevant_modified_crs <- st_transform(sf_relevant, crs = crs(proj.laea))
plot(sf_relevant_modified_crs)

# Come up with season cutoffs and trim data accordingly: use only non-breeding season data
# For Bjørnøya data:

bjo_nbs_df <- indiv_merged_df |> filter(colony == "Bjørnøya", !grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp)) |> dplyr::select(c(individ_id, timestamp, colony, lon, lat))
View(bjo_nbs_df)
bjo_sf_nbs_df <- st_as_sf(bjo_nbs_df, coords = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
length(unique(bjo_sf_nbs_df$individ_id))
plot(plastics)
plot(bjo_sf_nbs_df, cex = 0.4, col = "blue", pch = 16, add = T)
plot(crop(plastics, bjo_sf_nbs_df))
plot(bjo_sf_nbs_df, cex = 0.4, col = "blue", pch = 16, add = T)

# For Jan Mayen data: non-breeding season is defined here from April to September

jan_nbs_df <- indiv_merged_df |> filter(colony == "Jan Mayen", !grepl(c('-04-|-05-|-06-|-07-|-08-|-09-0') ,timestamp)) |> dplyr::select(c(individ_id, timestamp, colony, lon, lat))
jan_sf_nbs_df <- st_as_sf(jan_nbs_df, coords = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(plastics)
plot(jan_sf_nbs_df, cex = 0.4, col = "blue", pch = 16, add = T)
jan_nbs_cropped <- crop(plastics, jan_sf_nbs_df)
extent(jan_sf_nbs_df)
plot(jan_nbs_cropped, cex = 0.4, pch = 16)
jan_values_trial <- raster::extract(plastics, jan_sf_nbs_df, df = T)
length(unique(jan_sf_nbs_df$individ_id))

# To find: how to crop just vector point corresponding data from a raster

library(ggplot2)
na.omit(values(crop(plastics, bjo_sf_nbs_df)))

plastics_spdf <- as(plastics, "SpatialPixelsDataFrame")
plastics_spdf_df <- as.data.frame(plastics_spdf)
colnames(plastics_spdf_df) <- c("Value", "x", "y")

bjørnoya_overlap_gg <- ggplot() +  
  geom_tile(data = plastics_spdf_df, aes(x = x, y = y, fill = Value), alpha = 0.8) + 
  geom_sf(data = bjo_sf_nbs_df, alpha = 0.5, cex = 0.4) +
  scale_fill_viridis_c() +
  coord_sf() +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm")) +
  theme(axis.line = element_line("grey")) 
bjørnoya_overlap_gg

# Very ugly graph but shows points according to each of the 38 individuals

individ_bjørnoya_overlap_gg <- ggplot() +  
  geom_tile(data = plastics_spdf_df, aes(x = x, y = y, fill = Value), alpha = 0.8) + 
  geom_sf(data = bjo_sf_nbs_df, aes(col = individ_id), alpha = 0.5, cex = 0.4) +
  scale_fill_viridis_c() +
  coord_sf() +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm")) +
  theme(axis.line = element_line("grey")) 
individ_bjørnoya_overlap_gg

# Setting up a colour gradient for temporal ease of view - normal plotting: light blue denotes the
# start of the nbs (October, November), grey denotes intermediate months (December and January) and
# dark blue denotes the end of the nbs (February, March)

library(stringr)
plot(plastics)
plot(bjo_sf_nbs_df, 
     col = ifelse(str_detect(bjo_sf_nbs_df$timestamp, "-10-|-11-"), 5, 
                  ifelse(str_detect(bjo_sf_nbs_df$timestamp, "-02-|-03-"), 4, 8)), pch = 16, cex = 0.4, add = T)

# Recreating this in ggplot: change colour scaling - remove blue

temporal_bjo_overlap_gg <- ggplot() +  
  geom_tile(data = plastics_spdf_df, aes(x = x, y = y, fill = Value), alpha = 0.8) + 
  geom_sf(data = bjo_sf_nbs_df, aes( 
    col = ifelse(str_detect(timestamp, "-10-|-11-"), 1, 
                 ifelse(str_detect(timestamp, "-02-|-03-"), 5, 3))  
  ), alpha = 0.5, cex = 0.4) +
  scale_fill_viridis_c() +
  coord_sf() +
  theme_classic() +
  labs(fill = "Floating plastic debris value", col = "Time", 
       title = "Bjørnoya colony", 
       subtitle = "Northern fulmar locations during their non-breeding season overlayed on top of a global floating plastic debris distribution",
       tag = "A") +
  theme(legend.position="bottom", legend.key.width=unit(2, "cm"), axis.line = element_line("grey"),
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) 
temporal_bjo_overlap_gg

# Recreating the above graph for Jan Mayen

temporal_jan_overlap_gg <- ggplot() +  
  geom_tile(data = plastics_spdf_df, aes(x = x, y = y, fill = Value), alpha = 0.8) + 
  geom_sf(data = jan_sf_nbs_df, aes( 
    col = ifelse(str_detect(timestamp, "-10-|-11-"), 1, 
                 ifelse(str_detect(timestamp, "-02-|-03-"), 5, 3))  
  ), alpha = 0.5, cex = 0.4) +
  scale_fill_viridis_c() +
  coord_sf() +
  theme_classic() +
  labs(fill = "Floating plastic debris value", col = "Time", 
       title = "Jan Mayen colony", 
       subtitle = "Northern fulmar locations during their non-breeding season overlayed on top of a global floating plastic debris distribution",
       tag = "B") +
  theme(legend.position="bottom", legend.key.width=unit(2, "cm"), axis.line = element_line("grey"),
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) 
temporal_jan_overlap_gg

# Ignore for now- basic plot tis

jan_overlap_gg <- ggplot() +  
  geom_tile(data = plastics_spdf_df, aes(x = x, y = y, fill = Value), alpha = 0.8) + 
  geom_sf(data = jan_sf_nbs_df, alpha = 0.5, cex = 0.4) +
  scale_fill_viridis_c() +
  coord_sf() +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm")) +
  theme(axis.line = element_line("grey")) 
jan_overlap_gg

bjo_nbs_values <- raster::extract(x = plastics, y = bjo_sf_nbs_df, df = T)
jan_nbs_values <- raster::extract(x = plastics, y = jan_sf_nbs_df, df = T)

# Read phenology papers for stringent seasonal bounds; check petrels paper for what kind of analysis 
# they've done; what different analyses we could run, do plots like the ones in the fisheries paper where red and blue display temporal trends. 

# Checking out the adehabitat package in R

install.packages("adehabitatHR")
library(adehabitatHR)

trial_no <- relevant_new_data_1 |> dplyr::select(c('lon', 'lat'))
trial_no_matrix <- as.matrix(trial_no, ncol = 2)
trial_no_sp <- SpatialPoints(trial_no_matrix)
clu<- clusthr(trial_no_sp)
class(clu)
