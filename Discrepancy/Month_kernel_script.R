# This is a (relatively) clean script which calculates kernels on a month-wise basis, then sums months
# into one non-breeding season raster for each colony.

# Clearing environment
remove(list = ls())

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
library(spData) # for world vector
library(cowplot)
library(stringr)
library(RColorBrewer)
library(viridis)
library(viridisLite)


# Loading data----

datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/input_data/"
mylocs <- readRDS(paste0(datadir, "SEATRACK_FUGLA_20220307_v2.3_FA.rds"))
summary_info <- readRDS(paste0(datadir,"summaryTable.rds"))

# Exploring data
mylocs %>%
  dplyr::left_join(summary_info, join_by(ring)) %>%
  dplyr::mutate(year = year(timestamp)) %>%
  dplyr::group_by(colony,year) %>%
  dplyr::summarise(n_birds_tracked = n_distinct(ring))

# Creating an sf object----

indiv_merged_df <- merge(mylocs, summary_info, by = "ring") 
nbs_mylocs <- indiv_merged_df %>% 
  filter(!grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp)) %>%
  dplyr::mutate(year = year(timestamp))
names(nbs_mylocs)[names(nbs_mylocs) == "ring"] <- "individ_id"
nbs_mylocs_sf <- nbs_mylocs %>% 
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) 

# Assorting by months----

all_data <- nbs_mylocs
df <- all_data[!is.na(all_data$timestamp), ]
df <- df %>% dplyr::mutate(month = month(timestamp))
months <- sort(unique(df$month))

# Setting a directory for month-wise output rasters
dir_kernels <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/outputs"

# Defining land----

land <- as(world, "Spatial")

# Defining wgs84 projection
proj_wgs84 <- sp::CRS(sp::proj4string(land))

# In case, one wants to exactly replicate what the petrels' paper did
# setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Month/input_data/")
# land <- sf::read_sf(dsn = "baselayer", layer = "ne_10m_land")
# plot(land)
# proj_wgs84 <- crs(land)

# For loop starts here----

for(i in unique(df$colony)){ # First for loop begins
  sub <- df %>% filter(colony == i)

for(month_number in 1:length(months)){ # Second for loop begins
  tracks_wgs <- sub[sub$month == months[month_number],] 
  
sub <- subset(sub, lat < 90)  
png(filename = paste0(dir_kernels, "/distributions/", i, ".png"))  
plot(lat~lon, data = sub, type= "n", asp = 1, 
     xlim=c(min(all_data$lon) + 0.2,max(all_data$lon) - 0.2), 
     ylim=c(min(all_data$lat) + 0.5, max(all_data$lat) - 0.5), 
     main="", frame = T, xlab="", ylab="")
plot(land, col='lightgrey', add= T)
points(lat~lon, data = sub, pch = 16, cex = 0.5, col="blue")
points(lat~lon, data = tracks_wgs, pch = 16, cex = 0.5, col="red")
dev.off()

# Kernel density estimation----

if(nrow(tracks_wgs) > 4){
  
  if(min(tracks_wgs$lon) <= -179 ){ lon_min <- -180
  } else {lon_min <- floor(min(tracks_wgs$lon))-1 }
  
  if(max(tracks_wgs$lon) >= 179){ lon_max <- 180
  } else { lon_max <- ceiling(max(tracks_wgs$lon))+1 }
  
  if(min(tracks_wgs$lat) <= -89 ){ lat_min <- -90 
  } else { lat_min <- floor(min(tracks_wgs$lat))-1 }
  
  if(max(tracks_wgs$lat) >= 89){ lat_max <- 90
  } else { lat_max <- ceiling(max(tracks_wgs$lat))+1 }
  
  so.grid <- expand.grid( LON = seq(lon_min, lon_max, by=1), 
                          LAT = seq(lat_min, lat_max, by=1))
  
  sp::coordinates(so.grid) <- ~LON+LAT
  crs(so.grid) <- proj_wgs84
  
  # Setting a colony-centric crs
  mean_loc <- geosphere::geomean(cbind(tracks_wgs$lon,tracks_wgs$lat))
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
  sp::coordinates(null.grid) <- ~x+y
  sp::gridded(null.grid) <- TRUE

  # Converting tracks_wgs into a spatial points data frame
  sp::coordinates(tracks_wgs) <- ~lon+lat
  crs(tracks_wgs) <- proj_wgs84
  
  tracks <- sp::spTransform(tracks_wgs, CRS = DgProj)
  tracks$month <- factor(tracks@data$month)
  
  kudl <- adehabitatHR::kernelUD(tracks[,"month"], 
                                 grid = null.grid, h = 200000)  ## smoothing factor equals 200 km for GLS data
  
  vud <- adehabitatHR::getvolumeUD(kudl)
  
  fud <- vud[[1]]
  hr95 <- as.data.frame(fud)[,1]
  hr95 <- as.data.frame(fud)[,1]
  hr95 <- data.frame(hr95)
  coordinates(hr95) <- coordinates(fud)
  sp::gridded(hr95) <- TRUE
  
  kde_spixdf <- adehabitatHR::estUDm2spixdf(kudl)
  kern95 <- kde_spixdf
  
  stk_100 <- raster::stack(kern95)
  stk_95 <- raster::stack(hr95)
  
  sum_all_100 <- stk_100[[1]]
  sum_all_95 <- stk_95[[1]]
  
  sum_all_raw <- sum_all_100 * sum_all_95
  
  rast <- sum_all_raw/sum(raster::getValues(sum_all_raw))
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
  
  # Changing land's projection for masking
  mask_proj <- sp::spTransform(land, DgProj)
  mask_proj_pol <- as(mask_proj, "SpatialPolygons") 
  
  rast_mask_na <- raster::mask(cropped, mask_proj_pol, inverse = TRUE)
  rast_mask <- rast_mask_na
  rast_mask[is.na(rast_mask)] <- 0
  rast_mask_sum1 <- rast_mask/sum(raster::getValues(rast_mask))
  
  rast_mask_final <- raster::mask(rast_mask_sum1, mask_proj_pol, inverse = TRUE)
  rast_mask_final2 <- rast_mask_final 
  
  mask_wgs84 <- projectRaster(rast_mask_final2, crs = proj_wgs84, over = F)
  
  KDE_ref <- paste0(i, "_", months[month_number])
  raster::writeRaster(rast_mask_final, filename = paste0(datadir,"tifs/",KDE_ref,".tif"), 
                      format = "GTiff", overwrite = TRUE)
  raster::writeRaster(mask_wgs84, filename = paste0(dir_kernels,"/unique_tifs/",KDE_ref,".tif"), 
                      format = "GTiff", overwrite = TRUE)
  
  png(filename = paste0(dir_kernels,"/unique_distributions/",KDE_ref,".png"))
  plot(mask_wgs84, main = paste0("",i))
  plot(land, add = T, col = "#66000000")
  dev.off()

  
} # First if loop ends
} # Second for loop ends
} # First for loop ends



    
    
      
      
    
    
      
      
     
      
     
      
      
      
 
      
      
      
    
   
      
    
      
      
      
  
  


