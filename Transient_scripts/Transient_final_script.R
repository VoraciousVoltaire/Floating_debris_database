# Transient final script 

# KERNEL ESTIMATION----

# Loading relevant packages----
library(dplyr)
library(lubridate)
library(here)
library(ggplot2)
library(utils)
library(sf)
library(sp)
library(raster)
library(terra)
library(dplyr)
library(tidyverse)
library(adehabitatHR) 
library(amt)
library(rnaturalearth)
library(spData) # for world
library(cowplot)
library(stringr)
library(RColorBrewer)
library(viridis)
library(viridisLite)

# Loading og data
datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/input_data/fulmars_project_data/og_data/"
mylocs <- readRDS(paste0(datadir, "SEATRACK_FUGLA_20220307_v2.3_FA.rds"))
summary_info <- readRDS(paste0(datadir,"summaryTable.rds"))

# Exploring data
n_tracks <- mylocs %>%
  dplyr::left_join(summary_info, join_by(ring)) %>%
  dplyr::mutate(year = year(timestamp)) %>%
  dplyr::group_by(colony,year) %>%
  dplyr::summarise(n_birds_tracked = n_distinct(ring))
View(n_tracks)
total_n_tracks <- n_tracks %>% group_by(colony) %>% summarise(total = sum(n_birds_tracked))
View(total_n_tracks)

# Creating an sf object
indiv_merged_df <- merge(mylocs, summary_info, by = "ring")
nbs_mylocs <- indiv_merged_df %>% 
  filter(!grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp)) %>%
  dplyr::mutate(year = year(timestamp))
names(nbs_mylocs)[names(nbs_mylocs) == "ring"] <- "individ_id"
nbs_mylocs_sf <- nbs_mylocs %>% 
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) 

# Plotting nbs_mylocs
worldmap <- ggplot2::map_data('world')

ggplot() +
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group)) +
  geom_sf(data = nbs_mylocs_sf, aes(colour = colony), size = .1) +
  coord_sf(xlim = c(-90,90), ylim = c(10, 90)) +
  theme_light() +
  ggtitle("Northern Fulmar distribution in the non-breeding season")


# Defining land----
land <- as(world, "Spatial")

# Setting up a directory for outputs
dir_kernels <- paste0(getwd(), "/renditions_output/og_outputs")
proj_wgs84 <- sp::CRS(sp::proj4string(land))

# Kernel estimation prerequisites
# Setting a null grid----

if(min(nbs_mylocs$lon) <= -179 ){ lon_min <- -180
} else {lon_min <- floor(min(nbs_mylocs$lon))-1 }

if(max(nbs_mylocs$lon) >= 179){ lon_max <- 180
} else { lon_max <- ceiling(max(nbs_mylocs$lon))+1 }

if(min(nbs_mylocs$lat) <= -89 ){ lat_min <- -90 
} else { lat_min <- floor(min(nbs_mylocs$lat))-1 }

if(max(nbs_mylocs$lat) >= 89){ lat_max <- 90
} else { lat_max <- ceiling(max(nbs_mylocs$lat))+1 }

so.grid <- expand.grid(LON = seq(lon_min, lon_max, by=1), 
                       LAT = seq(lat_min, lat_max, by=1))

sp::coordinates(so.grid) <- ~LON+LAT
sp::proj4string(so.grid) <- sp::proj4string(land)

mean_loc <- geosphere::geomean(cbind(nbs_mylocs$lon, nbs_mylocs$lat))
DgProj <- sp::CRS(paste0("+proj=laea +lon_0=",mean_loc[1],
                         " +lat_0=",mean_loc[2])) 

so.grid.proj <- sp::spTransform(so.grid, CRS = DgProj)
coords <- so.grid.proj@coords

c <- min(coords[,1])-1000000   ## to check my min lon
d <- max(coords[,1])+1000000   ## to check my max lon

e <- min(coords[,2])-1000000   ## to check my min lat
f <- max(coords[,2])+1000000   ## to check my max lat

a <- seq(c, d, by=10000)
b <- seq(e, f, by=10000)
null.grid <- expand.grid(x = a,y = b)
sp::coordinates(null.grid) <- ~ x + y
sp::gridded(null.grid) <- TRUE

sp::coordinates(nbs_mylocs) <- ~ lon + lat
sp::proj4string(nbs_mylocs) <- sp::proj4string(land)

tracks <- sp::spTransform(nbs_mylocs, CRS = DgProj)
tracks$individ_id <- factor(tracks$individ_id)

# For loop----

cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

for(i in unique(nbs_mylocs$colony)){
  sub <- as.data.frame(nbs_mylocs) %>% filter(colony == i)
  sp::coordinates(sub) <- ~ lon + lat
  sp::proj4string(sub) <- sp::proj4string(land)
  tracks <- sp::spTransform(sub, CRS = DgProj)
  tracks$individ_id <- factor(tracks$individ_id)
  
  kudl <- kernelUD(tracks[,"individ_id"], grid = null.grid, h = 200000) ## smoothing factor equals 200 km for GLS data
  vud <- adehabitatHR::getvolumeUD(kudl)
  
  fud <- vector("list", length(unique(sub$individ_id)))
  hr95_cumulative <- data.frame(matrix(NA, nrow = 1)) 
  
  for(j in 1:length(unique(sub$individ_id))){
    fud <- vud[[j]]
    hr95 <- as.data.frame(fud)[,1]
    hr95 <- as.numeric(hr95 <= 95)
    hr95 <- data.frame(hr95)
    colnames(hr95) <- paste0(unique(sub$individ_id)[j])
    hr95_cumulative <- as.data.frame(cbind.fill(hr95_cumulative, hr95))
  }
  
  # Try this out with just the two colonies
 
  
  # Creating a cumulative spdf 
  hr95_mod <- hr95_cumulative[,-1]
  coordinates(hr95_mod) <- coordinates(fud)
  sp::gridded(hr95_mod) <- TRUE
  
  kern100 <- adehabitatHR::estUDm2spixdf(kudl)
  kern100_stk <- raster::stack(kern100)
  kern100_stk_mean <- terra::mean(kern100_stk, na.rm = T)
  
  volume95_stk <- raster::stack(hr95_mod)
  volume95_stk_mean <- terra::mean(volume95_stk, na.rm = T)
  
  sum_all_raw <- kern100_stk_mean * volume95_stk_mean
  
  rast <- sum_all_raw / sum(raster::getValues(sum_all_raw))
  rast[rast == 0] <- NA
  
  x.matrix <- is.na(as.matrix(rast))
  colNotNA <- which(colSums(x.matrix) != nrow(rast))
  rowNotNA <- which(rowSums(x.matrix) != ncol(rast))
  
  croppedExtent <- raster::extent(rast, 
                                  r1 = rowNotNA[1]-2, 
                                  r2 = rowNotNA[length(rowNotNA)]+2,
                                  c1 = colNotNA[1]-2, 
                                  c2 = colNotNA[length(colNotNA)]+2)
  
  cropped <- raster::crop(rast, croppedExtent)
  cropped[is.na(cropped)] <- 0
  
  # Masking land
  mask_proj <- sp::spTransform(land, DgProj) # changing projection 
  mask_proj_pol <- as(mask_proj, "SpatialPolygons") 
  
  ## set to NA cells that overlap mask (land)
  rast_mask_na <- raster::mask(cropped, mask_proj_pol, inverse = TRUE)
  rast_mask <- rast_mask_na
  rast_mask[is.na(rast_mask)] <- 0
  rast_mask_sum1 <- rast_mask/sum(raster::getValues(rast_mask))
  rast_mask[rast_mask == 0] <- NA
  rast_mask_final <- raster::mask(rast_mask_sum1, mask_proj_pol, inverse = TRUE)
  rast_mask_final2 <- rast_mask_final 
  
  #PLOT & SAVE ####
  mask_wgs84 <- raster::projectRaster(rast_mask_final2, crs = proj_wgs84, over = F)
  setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/renditions_output/og_outputs/")
  raster::writeRaster(mask_wgs84, filename = paste0("", i, ".tif"),
                      format = "GTiff", overwrite = T)
  
  ## Plot
  land <- as(world, "Spatial")
  land@bbox <- as.matrix(mask_wgs84@extent)
  png(filename = paste0("trial_", i, ".png"))
  plot(mask_wgs84)
  plot(land, add = T, col = "#66000000") # Next: Adjust scale over here to get apt images
  dev.off()
  
  View(hr95_cumulative)
  
}

# Multiplication rasters----

# Loading plastics raster
setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/input_data/fulmars_project_data/")
dir_1by1 <- ("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/renditions_output/og_outputs/")
plastics <- raster::raster("00_PlasticsRaster.tif")

r <- raster()
r

yelblus <- c(brewer.pal(n = 9, name = "YlGnBu"),"#00172e")
cols <- colorRampPalette(yelblus)(255)
colsviri <- cols[20:255]

colsinf <- rev(inferno(200))

## Rescale value sum to 1
plastics2 <- plastics
plastics2[is.na(plastics2)] <- 0
sum(raster::getValues(plastics2))
p_sum1    <- plastics2/sum(raster::getValues(plastics2)) # why, to make the values more generalizable?
p_sum1[is.na(plastics)] <- NA

RES <- res(plastics) # the resolution of the raster (in degrees)

R <- 6371007.2 # the Earth's authalic radius (in meters)- what is authalic???????????????????????????????????????????????????????????????
lat <- raster::yFromRow(plastics, 1:nrow(plastics)) # latitude of the centroid of each cell (in degrees, need to be converted in radians)
area <- (sin(pi/180*(lat + RES[2]/2)) - sin(pi/180*(lat - RES[2]/2))) * (RES[1] * pi/180) * R^2
r_area <- raster::setValues(plastics, rep(area, each = ncol(plastics))) # gives the area of each grid cell in meters 
plot(r_area, col = colsviri)

dir_demClasses <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/input_data/fulmars_project_data/og_data/Colony_rasters"
files <- list.files(dir_demClasses, full.names = TRUE, pattern = ".*\\.tif$"); files

dat <- data.frame()
result <- c()

nas <- as.data.frame(files)
nas$name <- NA
nas$vals <- NA
nas$nas <- NA
nas$files <- NULL

# Efficient for loop

for (i in 1:length(files)){
  
  a <- raster(files[i])
  name <- a@data@names[1]
  nas$name[i] <- name
  
  a[is.na(a)] <- 0 
  b <- median(raster::getValues(a)) 
  a_proj <- raster::projectRaster(a, plastics, method = "bilinear")
  a_proj2 <- a_proj * r_area / 100000000 # rescaling the values in each cell; but they haven't divided the values in the original plastics raster by their respective rasters areas or wasn't that necessary since the values are already normalized? Yeah, you have the proportions and you proportionally inflate the values according to the new projection system
  a_proj2[is.na(a_proj2)] <- 0 
  plot(a_proj2)
  c <- median(raster::getValues(a_proj2)) 
  
  # Find number of NAs for plastics where birds are----
  na_a_proj2 <- a_proj2
  na_a_proj2[na_a_proj2 > 0] <- 1
  na_a_proj2[is.na(na_a_proj2)] <- 0
  nas$vals[i] <- sum(na_a_proj2@data@values)
  
  na_p_sum1 <- p_sum1
  na_p_sum1[is.na(na_p_sum1)] <- 1
  na_p_sum1[na_p_sum1 != 1] <- 0
  na_over <- na_a_proj2 * na_p_sum1
  nas$nas[i] <- sum(na_over@data@values) 
  
  a_proj2[is.na(plastics)] <- NA
  
  ## Exporting results----
  raster_name_1 <- gsub(dir_demClasses, "", files[i])
  print(raster_name_1) 
  
  raster_name_2 <- paste0(dir_1by1, raster_name_1)
  
  raster::writeRaster(a_proj2, filename=raster_name_2, format="GTiff", overwrite=TRUE)
  
  over <- a_proj2 * p_sum1
  
  over_score <- over
  summary(over_score@data@values)
  over_score[is.na(over_score)] <- 0
  exposure_score <- round(sum(raster::getValues(over_score))*1000000, 4)
  
  png(paste0(dir_1by1, name,".png"), width= 1399,height= 455)
  par(mfrow=c(1,2))
  plot(a_proj2, main = paste0(name," distribution"),col= colsviri,legend = F)
  plot(over, main = paste0("Exposure score = ",exposure_score),
       col = colsinf,legend = F)
  dev.off()
  
}



# Note to myself: start writing about data collection like the loc_type column. Read and explain IRMA. Adjust scale in the png images to get appropriate zooming: land (scale = setextentshere)