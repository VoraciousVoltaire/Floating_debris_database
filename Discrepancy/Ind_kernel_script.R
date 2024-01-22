# This is a (relatively) clean script which calculates kernels on individual animal basis; calculates each individual's exposure score 
# and takes a mean of those values to arrive at one exposure score for each colony

# Clearing environment
rm(list = ls())

# Loading packages----

library(sf)
library(sp)
library(raster)
library(terra)
library(dplyr)
library(tidyverse)
library(adehabitatHR) # for calculating kernels and their volumes
library(amt)
library(spData) # for loading the spatial object 'world'

# Setting outputs directory
dir_kernels <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Ind/outputs"

# Loading data----

datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Ind/input_data"
mylocs <- readRDS(paste0(datadir, "/SEATRACK_FUGLA_20220307_v2.3_FA.rds"))
summary_info <- readRDS(paste0(datadir,"/summaryTable.rds"))

# Merging and curating for NBS (Non-Breeding Season)----

indiv_merged_df <- merge(mylocs, summary_info, by = "ring") 
nbs_mylocs <- indiv_merged_df %>% 
  filter(!grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp)) %>%
  dplyr::mutate(year = year(timestamp))
names(nbs_mylocs)[names(nbs_mylocs) == "ring"] <- "individ_id"

# Defining land----
land <- as(world, "Spatial")

# Defining wgs84 projection
proj_wgs84 <- sp::CRS(sp::proj4string(land))

# For loop for calculating individual animal kernels----

for(i in unique(nbs_mylocs$colony)){ # First for loop start
  sub <- as.data.frame(nbs_mylocs) %>% filter(colony == i) 
 
  # Setting a colony-centered crs 
  median_loc <- cbind(median(sub$lon), median(sub$lat))
  DgProj <- sp::CRS(paste0("+proj=laea +lon_0=",median_loc[1]," +lat_0=",median_loc[2]," +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m ", sep="")) 
 
  # Creating a colony-specific spatial points data frame
  sp::coordinates(sub) <- ~ lon + lat
  sp::proj4string(sub) <- sp::proj4string(land)
  
  # Setting a null grid for kernel estimation 
  if(nrow(sub) > 4){
    
    if(min(sub$lon) <= -179 ){ lon_min <- -180
    } else {lon_min <- floor(min(sub$lon))-1 }
    
    if(max(sub$lon) >= 179){ lon_max <- 180
    } else { lon_max <- ceiling(max(sub$lon))+1 }
    
    if(min(sub$lat) <= -89 ){ lat_min <- -90 
    } else { lat_min <- floor(min(sub$lat))-1 }
    
    if(max(sub$lat) >= 89){ lat_max <- 90
    } else { lat_max <- ceiling(max(sub$lat))+1 }
    
    so.grid <- expand.grid(LON = seq(lon_min, lon_max, by=1), 
                           LAT = seq(lat_min, lat_max, by=1))
    
    sp::coordinates(so.grid) <- ~LON+LAT
    sp::proj4string(so.grid) <- sp::proj4string(land)
    
    so.grid.proj <- sp::spTransform(so.grid, CRS = DgProj) # transforming null grid to new crs
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
  
  tracks <- sp::spTransform(sub, CRS = DgProj)
  tracks$individ_id <- factor(tracks$individ_id)
  
# Kernel estimation
  kudl <- kernelUD(tracks[,"individ_id"], grid = null.grid, h = 200000) ## smoothing factor equals 200 km for GLS data
  
# Changing from estUDm class to spatial pixels data frame
  kern100 <- adehabitatHR::estUDm2spixdf(kudl)
# Creating a raster stack so that each layer corresponds to a raster from a separate individual
  kern100_stk <- raster::stack(kern100)
  
  vud <- adehabitatHR::getvolumeUD(kudl) # calculating home range percent estimates from kernels
  
  fud <- vector("list", length(unique(sub$individ_id)))
  
  for(j in 1:length(unique(sub$individ_id))){ # Second for loop starts
    fud <- vud[[j]]
    hr95 <- as.data.frame(fud)[,1]
    hr95 <- as.numeric(hr95 <= 95) # now the data's binary, 1 meaning the pixel belongs to the individual's home range
    hr95 <- data.frame(hr95)
  
  # Creating a spatial pixels data frame
  hr95_mod <- hr95
  coordinates(hr95_mod) <- coordinates(fud)
  sp::gridded(hr95_mod) <- TRUE
  hr95_raster <- raster(hr95_mod)


  rast <- kern100_stk[[j]] * hr95_raster
  rast[is.na(rast)] <- 0
    
    # Proportionate scaling of home range rasters
    rast <- rast/sum(raster::getValues(rast))
    
    # Cropping extent
    rast[rast == 0] <- NA
    x.matrix <- is.na(as.matrix(rast))
    colNotNA <- which(colSums(x.matrix) != nrow(rast))
    rowNotNA <- which(rowSums(x.matrix) != ncol(rast))
    croppedExtent <- raster::extent(rast, 
                                    rowNotNA[1] - 2, 
                                    rowNotNA[length(rowNotNA)] + 2,
                                    colNotNA[1] - 2, 
                                    colNotNA[length(colNotNA)] + 2)
    cropped <- raster::crop(rast, croppedExtent)
    cropped[is.na(rast)] <- 0
    
    # Masking land
    mask_proj <- sp::spTransform(land, DgProj) # changing projection 
    mask_proj_pol <- as(mask_proj, "SpatialPolygons") 
    
    ## set to NA cells that overlap mask (land)
    rast_mask_na <- raster::mask(cropped, mask_proj_pol, inverse = TRUE)
    rast_mask <- rast_mask_na
    rast_mask[is.na(rast_mask)] <- 0
    rast_mask_sum1 <- rast_mask/sum(raster::getValues(rast_mask))
    rast_mask[rast_mask == 0] <- NA
    rast_mask_final <- raster::mask(rast_mask_sum1, mask_proj_pol, inverse = TRUE) # this is to be used as the input whilst doing raster multiplication
    
    
    #PLOT & SAVE ####
    mask_wgs84 <- raster::projectRaster(rast_mask_final,  crs = proj_wgs84, over = F)
    setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Ind/outputs/unique_tifs/")
    raster::writeRaster(mask_wgs84, filename = paste0(i,"_", j ,".tif"),
                        format = "GTiff", overwrite = T)      
    setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Ind/input_data/tifs/")
    raster::writeRaster(rast_mask_final, filename = paste0(i,"_", j ,".tif"),
                        format = "GTiff", overwrite = T)    
                     
    ## Plot
    png(filename = paste0(dir_kernels, "/unique_distributions/", i, "_", j, ".png"))
    plot(mask_wgs84)
    plot(land, add = T, col = "#66000000")
    dev.off()
    
  
    } # Second for loop ends
  
  } # first if loop end

} # first for loop end




