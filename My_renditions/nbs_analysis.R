sessionInfo()

# Loading relevant packages
library(sf)
library(sp)
library(raster)
library(terra)
library(dplyr)
library(tidyverse)
library(adehabitatHR) 
library(amt)
library(rnaturalearth)


# Working with non-breeding season's data----

setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/input_data/fulmars_project_data/")
new_data_1 <- readRDS("test_2colonies.rds")
new_data_2 <- readRDS("test_2colonies_individ_info.rds")
indiv_merged_df <- merge(new_data_1, new_data_2, by = "individ_id")
relevant_new_data_1 <- dplyr::select(indiv_merged_df, individ_id, timestamp, lon, lat, loc_type, colony)
df <- st_as_sf(relevant_new_data_1, coords = c('lon','lat'), crs = 4326)
relevant_new_data_2 <- relevant_new_data_1 %>% filter(!grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp))

# Condensed script----
# setting an overarching for loop for looping all colonies; using functions from adehabitatHR package
setwd('/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/')
library(adehabitatHR)
for(i in unique(relevant_new_data_1$colony)){
  sub <- relevant_new_data_1 %>% filter(colony == i, !grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp))
  sub_sf <- st_as_sf(sub, coords = c("lon", "lat"), crs = 4326)
  sub_coords <- sub[,c("lon", "lat")]
  sub_spdf <- SpatialPointsDataFrame(coords = sub_coords, data = sub, proj4string = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))
  sub_crs <-CRS(paste("+proj=laea +lat_0=",median(sub$lat)," +lon_0=",median(sub$lon)," +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km ", sep = ""))
  sub_spdf <- spTransform(sub_spdf, CRS = sub_crs)
  kernel <- kernelUD(sub_spdf[,1])
  png(paste0('renditions_output/loop_outputs/kernel_image_',i,'.png'), width = 1379, height = 750)
  plot(unlist(image(kernel)))
  ud <- getverticeshr(kernel, percent = 95)
  print(ud)
  ud_sf <- st_as_sf(ud)
  png(paste0('renditions_output/loop_outputs/UD_image_',i,'.png'), width = 1379, height = 750)
  plot(st_geometry(ud_sf[1,]))
  png(paste0('renditions_output/loop_outputs/coloured_UD_image',i,'.png'), width = 1379, height = 750)
  ud@data$id <- as.factor(ud@data$id)
  plot(ud, col = ud@data$id)
  dev.off()
}

# setting an overarching for loop for looping all colonies; using functions from amt package
# coupled_countries_all <- ne_countries(scale = 50, country= c("Norway", "Sweden", "Finland", "Russia", "Greenland", "United Kingdom", "Iceland", "Ireland", "Faroe Islands", "Denmark", "Netherlands", "Belgium", "France", "Canada","United States"), returnclass = 'sf')
# x <- seq(-5000,5000, by=1.) # resolution is the pixel size you desire 
# y <- seq(-5000, 5000, by=1.)
# xy <- expand.grid(x=x,y=y)
# coordinates(xy) <- ~x+y
# gridded(xy) <- TRUE
for(i in unique(relevant_new_data_1$colony)){
  sub <- relevant_new_data_1 %>% filter(colony == i, !grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp))
  # medlat <- median(sub$lat)
  # medlon <- median(sub$lon)
  # sub <- st_as_sf(sub, coords = c('lon', 'lat'), crs = 4326) 
  # sub <- st_transform(sub, crs = CRS(paste("+proj=laea +lat_0=",medlat," +lon_0=",medlon," +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km ", sep = "")))
  # sub <- sub %>% mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2])
  # sub <- st_set_geometry(sub, NULL)
  # sub_crs <-CRS(paste("+proj=laea +lat_0=",median(sub$lat)," +lon_0=",median(sub$lon)," +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km ", sep = ""))
  dat.track <- make_track(sub, lon, lat, timestamp, crs = 4326, all_cols = T)
  dat.mcp <- hr_mcp(dat.track, levels = c(0.5,0.95,1))
  ggplot() +
    geom_sf(data = ne_countries(scale = 50)) +
    geom_point(data = sub, aes(lon, lat), alpha = 0.05, size = 1) +
    geom_sf(data = dat.mcp$mcp, aes(colour = factor(level)), fill = 'transparent', size = 0.75) +
    scale_colour_viridis_d(direction = -1) +
    theme_bw() +
    coord_sf(xlim = c(min(sub$lon) - 5, max(sub$lon) + 5), ylim = c(min(sub$lat) - 5, max(sub$lat) + 5)) # limits can be changed later 
  ggsave(paste0('renditions_output/loop_outputs/coloured_mcp_',i,'.png'))
  trast <- make_trast(dat.track)
  h_pi <- hr_kde_pi(dat.track)
  dat.kde.pi <- hr_kde(dat.track, trast = trast, h = h_pi, levels = c(0.5,0.95))
  # plot(dat.kde.pi, col = c("red", "blue"))
  # png(paste0('renditions_output/loop_outputs/amt_kde_',i,'.png')) 
  # Warning message:
  # In KernSmooth::bkde2D(as.matrix(x[, c("x_", "y_")]), bandwidth = h,  :
  #                         Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'
  kde.hpi.contours <- hr_isopleths(dat.kde.pi)
  ggplot() +
    geom_sf(data = ne_countries(scale = 50)) +
    geom_path(data = data.frame(sub), aes(lon, lat, group = individ_id), alpha = 0.25, size = 0.3) +
    geom_sf(data = kde.hpi.contours, aes(colour = factor(level)), fill = 'transparent', size = 0.75) +
    scale_fill_viridis_c() +
    coord_sf(xlim = c(min(sub$lon) - 5, max(sub$lon) + 5), ylim = c(min(sub$lat) - 5, max(sub$lat) + 5)) # limits can be changed later
  ggsave(paste0('renditions_output/loop_outputs/amt_all_ind_kde_',i,'.png'))
  dat.id.kde.hpi <- dat.track %>%
    split(.$individ_id) %>%
    map(~hr_kde(.x,
                trast = make_trast(.x),
                h = hr_kde_pi(.x),
                levels = c(0.5, 0.95))) %>%
    map(hr_isopleths) %>% 
    do.call(rbind, .)
  dat.id.kde.hpi <- dat.id.kde.hpi %>%
    mutate(individ_id = rownames(.), .before = level)
  ggplot() +
    geom_sf(data = ne_countries(scale = 50)) +
    geom_path(data = sub, aes(lon, lat, group = individ_id), alpha = 0.25, size = 0.3) +
    geom_sf(data = dat.id.kde.hpi, aes(color = factor(level)), fill = 'transparent', size = 0.75) +
    scale_fill_viridis_c() +
    coord_sf(xlim = c(min(sub$lon) - 10, max(sub$lon) + 10), ylim = c(min(sub$lat) - 10, max(sub$lat) + 10)) +
    facet_wrap(~individ_id)
  ggsave(paste0('renditions_output/loop_outputs/amt_unique_ind_kde_',i,'.png'))
}

sf_relevant <- st_as_sf(relevant_new_data_1[,c(3,4,6)], coords = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Come up with season cutoffs and trim data accordingly: use only non-breeding season data
# For Bjørnøya data:

bjo_nbs_df <- indiv_merged_df |> filter(colony == "Bjørnøya", !grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp)) |> dplyr::select(c(individ_id, timestamp, colony, lon, lat))
View(bjo_nbs_df)
bjo_nbs_df[,1] <- as.factor(bjo_nbs_df[,1])
bjo_sf_nbs_df <- st_as_sf(bjo_nbs_df, coords = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
length(unique(bjo_sf_nbs_df$individ_id)) # 38 unique individuals
bjo_values <- raster::extract(plastics, bjo_sf_nbs_df, df = T)

# For Jan Mayen data: non-breeding season is defined here from April to September

jan_nbs_df <- indiv_merged_df |> filter(colony == "Jan Mayen", !grepl(c('-04-|-05-|-06-|-07-|-08-|-09-0') ,timestamp)) |> dplyr::select(c(individ_id, timestamp, colony, lon, lat))
jan_nbs_df[,1] <- as.factor(jan_nbs_df[,1])
jan_sf_nbs_df <- st_as_sf(jan_nbs_df, coords = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
length(unique(jan_sf_nbs_df$individ_id)) # 48 unique individuals

# Trying adehabitatHR for home range estimation: kernel density estimation was used in the 2023
# petrels paper so might be interesting to dig in more

# Checking out the adehabitat package----

install.packages("adehabitatHR")
library(adehabitatHR)

bjo_ade <- bjo_nbs_df |> dplyr::select(c('lon', 'lat'))
bjo_ade_matrix <- as.matrix(bjo_ade, ncol = 2)
bjo_ade_sp <- SpatialPoints(bjo_ade_matrix)
cp <- mcp(bjo_ade_sp)
as.data.frame(cp)

# Estimating home range for each individual using minimum convex polygon (MCP)----
# Creating a SpatialPointsDataFrame 
bjo_coords <- bjo_nbs_df[,c("lon","lat")]
jan_coords <- jan_nbs_df[,c("lon","lat")]

# Assigning bjo's geomedian-based crs to bjo_spdf

bjo_spdf <- SpatialPointsDataFrame(coords = bjo_coords, data = bjo_nbs_df, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
bjo_crs <-CRS("+proj=laea +lat_0=73.42981 +lon_0=20.8018574294084 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km ")
bjo_spdf_1 <- spTransform(bjo_spdf, bjo_crs)
bjo_cp <- mcp(bjo_spdf_1[,1])
as.data.frame(bjo_cp)
plot(bjo_cp)
plot(bjo_spdf_1, col=as.data.frame(bjo_spdf_1)[,1], add = T) 
bjo_cp
diff_lev <- mcp.area(bjo_spdf_1[,1]) # home range across different levels
kernel <- kernelUD(bjo_spdf_1[,1]) # finally, kernel densities!
image(kernel)
head(kernel)
image(kernel[[1]])
ud <- getverticeshr(kernel, percent = 90)
ud_1 <- st_as_sf(ud)
plot(st_geometry(ud_1[1,]))
# messed up for now, will find out how to segregate each of them 38 polygons

# Okay for subsetting unique individuals use- set a loop to get kernel densities for each individual
# First create a list of all unique individuals
bjo_unique_ind <- unique(bjo_spdf$individ_id)
bjo_spdf_df <- as.data.frame(bjo_spdf)
grouped_bjo <- bjo_spdf_df |> group_by(individ_id) |> summarize(n = n())
View(grouped_bjo)

output_df <- list()
bjo_spdf_df <- as.data.frame(bjo_spdf)
for(i in 1:length(bjo_unique_ind)){
  output_df[[i]] <- bjo_spdf_df |> filter(individ_id == bjo_unique_ind[i])
}
output_df

# I'd first have to convert output_df into SpatialPointsDataFrames
spatial_output_df <- list()
bjo_spatial_coords <- list()
for(i in 1:length(bjo_unique_ind)){
  bjo_spatial_coords[[i]] <- output_df[[i]][,c("lon","lat")]
  spatial_output_df[[i]] <- SpatialPointsDataFrame(coords = bjo_spatial_coords[[i]], data = output_df[[i]], proj4string = bjo_crs)
}

# Setting up another for loop for calculating 95% UDs
x <- seq(-75, 175, by=1.) # resolution is the pixel size you desire 
y <- seq(0, 100, by=1.)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE
class(xy)

bjo_kernel <- vector("list", 38)
bjo_ud <- list("list", 38)
for(i in 1:length(bjo_unique_ind)){
  bjo_kernel[[i]] <- kernelUD(spatial_output_df[[i]][,3], grid = xy)
}

  for (i in 1:length(bjo_unique_ind)){
    bjo_ud[[i]] <- getverticeshr(bjo_kernel[[i]], percent = 95)
  }

NOS-4181262
subset_1 <- subset(bjo_spdf, individ_id == "NOS-4181262")
unique(bjo_spdf$individ_id)
head(subset_1)
subset_1_kernel <- kernelUD(subset_1[,3])
image(subset_1_kernel)
subset_1_ud <- getverticeshr(subset_1_kernel, percent = 90)
plot(subset_1_ud, add = T)

# # trying out stuff from documentation
# volume <- getvolumeUD(subset_1_kernel)
# levels <- c(50, 75, 90)
# listt <- vector(mode = "list", length = 2)
# listt[[1]] <- as.image.SpatialGridDataFrame(volume[[1]])
# contour(list[[1]], levels = levels)

head(subset_1_ud)
library(sf)
fort <- fortify(subset_1_ud)
head(fort)
ggplot(fort, aes(x = long, y = lat)) + geom_polygon(alpha = 0.4) + coord_equal()

# ggtiles as the plastic raster
library(ggplot2)
plastics_spdf <- as(plastics, "SpatialPixelsDataFrame")
plastics_spdf_df <- as.data.frame(plastics_spdf)
colnames(plastics_spdf_df) <- c("Value", "x", "y")

bjo_kernel_overlap <- ggplot(fort, aes(x = long, y = lat)) +  
  coord_equal() +
  geom_tile(data = plastics_spdf_df, aes(x = x, y = y, fill = Value), alpha = 0.8) + 
  scale_fill_viridis_c() +
  geom_polygon(alpha = 0.6) +
  # theme_classic() +
  theme(legend.position="bottom", legend.key.width=unit(2, "cm"), axis.line = element_line("grey"),
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  labs(fill = "Floating plastic debris value",  
       title = "Bjørnoya colony", 
       subtitle = "95% utilization distributions of northern fulmars from Bjørnoya colony during their non-breeding season",
       tag = "A")
bjo_kernel_overlap

# Trials
install.packages("amt")
library(amt)
library(units) # to show spatial units in ggplots
library(MetBrewer) # to display colour palettes
install.packages("rnaturalearth")
library(rnaturalearth)
install.packages("aniMotum")

# Wrangling and preparing data for KDE
dat.track <- make_track(bjo_nbs_df, lon, lat, crs = "+proj=laea +lat_0=73.42981 +lon_0=20.8018574294084 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km")
# Template raster for KDE to map on to
trast <- make_trast(dat.track, res = 10) 
trast
# Calculating kde for individual ids
dat.kde.ref <- hr_kde(dat.track, trast = trast, h = hr_kde_ref(dat.track), levels = c(0.5,0.95))
dat.kde.ref
plot(dat.kde.ref, col = c("red", "blue"))
kde.href.contours <- hr_isopleths(dat.kde.ref)

# ggplot() + geom_tile(data = raster::as.data.frame(dat.kde.ref$ud)) +
#   scale_fill_viridis_c(option = "inferno") +
#   theme_bw()



NOS-4181268
subset_2 <- subset(bjo_spdf, individ_id == "NOS-4181268")
subset_2_kernel <- kernelUD(subset_2[,3])
image(subset_2_kernel)
subset_2_ud <- getverticeshr(subset_2_kernel, percent = 90)
plot(subset_2_ud, add = T)

NOS-5109523
subset_38 <- subset(bjo_spdf, individ_id == "NOS-5109523")
subset_38_kernel <- kernelUD(subset_38[,3])
image(subset_38_kernel)
subset_38_ud <- getverticeshr(subset_38_kernel, percent = 95)
plot(subset_38_ud, add = T)

bjo_unique_ind


jan_coords <- jan_nbs_df[,c("lon","lat")]
jan_spdf <- SpatialPointsDataFrame(coords = jan_coords, data = jan_nbs_df)
jan_cp <- mcp(jan_spdf[,1])
as.data.frame(jan_cp)
plot(jan_cp)

clu <- clusthr(bjo_ade_sp) # takes too long to execute
class(clu)
ade_df <- bjo_nbs_df[,-3] 
ade_spdf <- SpatialPointsDataFrame(coords = c(ade_df$lon, ade_df$lat), data = ade_df)
cp <- mcp(ade_df[,3])
?SpatialPointsDataFrame

# Setting a new crs based on median positions of geolocations---- 
# no utility as of now

# For Bjørnoya 
bjo_medlat = median(bjo_nbs_df$lat)
bjo_medlon = median(bjo_nbs_df$lon)
bjo_proj.laea = paste("+proj=laea +lat_0=",bjo_medlat, " +lon_0=",bjo_medlon," +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km ", sep="")
bjo_mod <- st_transform(bjo_sf_nbs_df, crs = crs(bjo_proj.laea))
plot(bjo_mod$geometry)

# For Jan Mayen
jan_medlat = median(jan_nbs_df$lat)
jan_medlon = median(jan_nbs_df$lon)
jan_proj.laea = paste("+proj=laea +lat_0=",jan_medlat, " +lon_0=",jan_medlon," +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km ", sep="")
jan_mod <- st_transform(jan_sf_nbs_df, crs = crs(jan_proj.laea))
plot(jan_mod$geometry)

# Plotting graphs for clarity 

# Bjørnoya
plot(plastics)
plot(bjo_sf_nbs_df, cex = 0.4, col = "blue", pch = 16, add = T)
plot(crop(plastics, bjo_sf_nbs_df))
plot(bjo_sf_nbs_df, cex = 0.4, col = "blue", pch = 16, add = T)

# Jan Mayen
plot(plastics)
plot(jan_sf_nbs_df, cex = 0.4, col = "blue", pch = 16, add = T)
jan_nbs_cropped <- crop(plastics, jan_sf_nbs_df)
extent(jan_sf_nbs_df)
plot(jan_nbs_cropped, cex = 0.4, pch = 16)

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


