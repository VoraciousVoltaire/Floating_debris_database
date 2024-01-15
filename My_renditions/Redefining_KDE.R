# Redefining kernel density estimations----

# Clearing environment----
rm(list = ls())

# Loading relevant packages----
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

# Loading data----
setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/input_data/fulmars_project_data/")
new_data_1 <- readRDS("test_2colonies.rds")
new_data_2 <- readRDS("test_2colonies_individ_info.rds")
indiv_merged_df <- merge(new_data_1, new_data_2, by = "individ_id")
relevant_new_data_1 <- dplyr::select(indiv_merged_df, individ_id, timestamp, lon, lat, loc_type, colony)
df <- st_as_sf(relevant_new_data_1, coords = c('lon','lat'), crs = 4326)
relevant_new_data_2 <- relevant_new_data_1 %>% filter(!grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp))

# Defining land----
land <- as(world, "Spatial")

# Setting up a directory and object for feasibility----
dir_kernels <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/renditions_output"
proj_wgs84 <- sp::CRS(sp::proj4string(land))

# Setting a null grid----

if(min(relevant_new_data_2$lon) <= -179 ){ lon_min <- -180
} else {lon_min <- floor(min(relevant_new_data_2$lon))-1 }

if(max(relevant_new_data_2$lon) >= 179){ lon_max <- 180
} else { lon_max <- ceiling(max(relevant_new_data_2$lon))+1 }

if(min(relevant_new_data_2$lat) <= -89 ){ lat_min <- -90 
} else { lat_min <- floor(min(relevant_new_data_2$lat))-1 }

if(max(relevant_new_data_2$lat) >= 89){ lat_max <- 90
} else { lat_max <- ceiling(max(relevant_new_data_2$lat))+1 }

so.grid <- expand.grid(LON = seq(lon_min, lon_max, by=1), 
                       LAT = seq(lat_min, lat_max, by=1))

sp::coordinates(so.grid) <- ~LON+LAT
sp::proj4string(so.grid) <- sp::proj4string(land)

mean_loc <- geosphere::geomean(cbind(relevant_new_data_2$lon, relevant_new_data_2$lat))
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

sp::coordinates(relevant_new_data_2) <- ~ lon + lat
sp::proj4string(relevant_new_data_2) <- sp::proj4string(land)

tracks <- sp::spTransform(relevant_new_data_2, CRS = DgProj)
tracks$individ_id <- factor(tracks$individ_id)

# For loop for calculating kde----

# Working loop!!!!!!!!!!!!!!!!!
setwd(dir_kernels)


for(i in unique(relevant_new_data_2$colony)){
  sub <- as.data.frame(relevant_new_data_2) %>% filter(colony == i)
  sp::coordinates(sub) <- ~ lon + lat
  sp::proj4string(sub) <- sp::proj4string(land)
  tracks <- sp::spTransform(sub, CRS = DgProj)
  tracks$individ_id <- factor(tracks$individ_id)
  
  kudl <- kernelUD(tracks[,"individ_id"], grid = null.grid, h = 200000) ## smoothing factor equals 200 km for GLS data
  vud <- adehabitatHR::getvolumeUD(kudl)
  
  fud <- vector("list", length(unique(sub$individ_id)))
  hr95_cumulative <- data.frame(matrix(NA, nrow = 932694))
  
  for(j in 1:length(unique(sub$individ_id))){
    fud <- vud[[j]]
    hr95 <- as.data.frame(fud)[,1]
    hr95 <- as.numeric(hr95 <= 95)
    hr95 <- data.frame(hr95)
    colnames(hr95) <- paste0(unique(sub$individ_id)[j])
    hr95_cumulative <- cbind(hr95_cumulative, hr95)
  }
  
    hr95_mod <- hr95_cumulative[,-1]
    coordinates(hr95_mod) <- coordinates(fud)
    sp::gridded(hr95_mod) <- TRUE
  
  
  #} # temporary break 
  
  # note
  # sample <- as.vector(unlist(vud))
  # range(vud[[1]]$n)
  # #note
  
  
  # Additions now
  
  setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/renditions_output/")
  
  # for(j in 1:length(kudl)){
  #   par(mfrow = c(2,1))
  #   xyz <- as.image.SpatialGridDataFrame(kudl[[j]])
  #   xyzv <- as.image.SpatialGridDataFrame(vud[[j]])
  #   # plot(unlist(image(kudl[[j]])))
  #  
  #   image(kudl[[j]])
  #   title("Output of kernelUD")
  #   contour(xyz, add = T)
  #  
  #   image(vud[[j]])
  #   title("Output of getvolumeUD")
  #   contour(xyzv, add = T)
  # }
  
  # Additions end

  # fud <- vector("list", length(unique(sub$individ_id)))
  # hr95 <- vector("list", length(unique(sub$individ_id)))
  # for(p in 1:length(unique(sub$individ_id))){
  #   class(vud[[p]])
  #   hr95 <- as.data.frame(fud[[p]][,1])
  # }
  # View(as.data.frame(estUDm2spixdf(vud)[,2]))
  
  
  #  # didn't understand why; shall figure out
  # hr95 <- as.data.frame(fud)[,1]
  # hr95 <- as.numeric(hr95 <= 95)
  # hr95 <- data.frame(hr95)
  # coordinates(hr95) <- coordinates(fud) 
  # sp::gridded(hr95) <- TRUE
  # 
  # # Additions start here
  # 
  # image(hr95)
  
  # I still have to figure out how to modify my new hr95 to get a 95% home range for each individual. Shouldn't be too hard, just have to find the right argument. 
  # And it's crucial to change it from 100 to 95 so that externalities are dealt with in a better manner
  
  # Additions end here
  
  volume100 <- adehabitatHR:: estUDm2spixdf(vud)
  kern100 <- adehabitatHR::estUDm2spixdf(kudl)
    
  kern100_stk <- raster::stack(kern100)
  kern100_stk_mean <- terra::mean(kern100_stk, na.rm = T) # for some reason, median wasn't working
  
  volume95_stk <- raster::stack(hr95_mod)
  volume95_stk_mean <- terra::mean(volume95_stk, na.rm = T)
  
  volume100_stk <- raster::stack(volume100)
  volume100_stk_mean <- terra::mean(volume100_stk, na.rm = T)
  
  sum_all_raw <- kern100_stk_mean * volume95_stk_mean
  
  # sum_all_raw <- sum_all_100[[k]] * sum_all_95[[k]]
  # s <- as.data.frame(as.vector(values(sum_all_raw)))
  # sum(s)
  
  
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
  setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/renditions_output/revised_script_12_1/")
  raster::writeRaster(mask_wgs84, filename = paste0("", i, ".tif"),
                      format = "GTiff", overwrite = T)
  
  ## Plot
  png(filename = paste0("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/renditions_output/revised_script_12_1/trial_", i, ".png"))
  plot(mask_wgs84)
  plot(land, add = T, col = "#66000000")
  dev.off()
}









  

# Investigating if smoothing factor can be manipulated in adhabitatHR
# Since I'm getting an error which says that the grid is too small to allow the estimation
# of home-range, you should rerun kernelUD with a larger extent parameter; I'll try 
# setting up a null grid as the petrels paper did

# # Setting a null grid (just copy pasted the previous code)----

if(min(relevant_new_data_2$lon) <= -179 ){ lon_min <- -180
} else {lon_min <- floor(min(relevant_new_data_2$lon))-1 }

if(max(relevant_new_data_2$lon) >= 179){ lon_max <- 180
} else { lon_max <- ceiling(max(relevant_new_data_2$lon))+1 }

if(min(relevant_new_data_2$lat) <= -89 ){ lat_min <- -90 
} else { lat_min <- floor(min(relevant_new_data_2$lat))-1 }

if(max(relevant_new_data_2$lat) >= 89){ lat_max <- 90
} else { lat_max <- ceiling(max(relevant_new_data_2$lat))+1 }

so.grid <- expand.grid(LON = seq(lon_min, lon_max, by=1), 
                       LAT = seq(lat_min, lat_max, by=1))

sp::coordinates(so.grid) <- ~LON+LAT
sp::proj4string(so.grid) <- sp::proj4string(land)

mean_loc <- geosphere::geomean(cbind(relevant_new_data_2$lon, relevant_new_data_2$lat))
DgProj <- sp::CRS(paste0("+proj=laea +lon_0=",mean_loc[1],
                         " +lat_0=",mean_loc[2])) 

so.grid.proj <- sp::spTransform(so.grid, CRS=DgProj)
coords <- so.grid.proj@coords

c <- min(coords[,1])-1000000   ## to check my min lon
d <- max(coords[,1])+1000000   ## to check my max lon

e <- min(coords[,2])-1000000   ## to check my min lat
f <- max(coords[,2])+1000000   ## to check my max lat

a <- seq(c, d, by=10000)
b <- seq(e, f, by=10000)
null.grid <- expand.grid(x=a,y=b)
sp::coordinates(null.grid) <- ~x+y
sp::gridded(null.grid) <- TRUE

# Changing CRS----
tracks <- sp::spTransform(relevant_new_data_2, CRS=DgProj)
tracks$individ_id <- factor(tracks@data$individ_id)

#Trying out plastics spdf as grid----
plastics_sp <- as(plastics, "SpatialPixels")
plastics_spdf_df <- as.data.frame(plastics_spdf)
colnames(plastics_spdf_df) <- c("Value", "x", "y")
plastics_2 <- plastics

setwd('/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/')
for(i in unique(relevant_new_data_1$colony)){
  sub <- relevant_new_data_1 %>% filter(colony == i, !grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp))
  sub_sf <- st_as_sf(sub, coords = c("lon", "lat"), crs = 4326)
  sub_coords <- sub[,c("lon", "lat")]
  sub_crs <-CRS(paste("+proj=laea +lat_0=",median(sub$lat)," +lon_0=",median(sub$lon)," +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km ", sep = ""))
  sub_spdf <- SpatialPointsDataFrame(coords = sub_coords, data = sub, proj4string = sub_crs)
  kernel <- kernelUD(sub_spdf[,1])
  
# EXPERIMENTATION STARTS HERE
  
  # # sub_SpatialPoints <- SpatialPoints(sub_spdf[3:4])
  # kernel <- kernelUD(sub_spdf[,1], grid = plastics_sp) # I'm making this happen, to have a kernel displayed over the plastics grid to make things more feasible?
  # 
  # 
  # sub_raster <- raster(sub_spdf)
  # # Reprojecting sub_raster
  # 
  # # The whole deal about plastics_2 and transformed plastics----
  # plastics_2 <- plastics
  # transformed_plastics <- projectRaster(plastics_2, crs = sub_crs)
  # 
  # # Prep for ggplot
  # trasnformed_plastics_spdf <- as(transformed_plastics, 'SpatialPixelsDataFrame')
  # transformed_plastics_spdf_df <- as.data.frame(transformed_plastics_spdf)
  # colnames(transformed_plastics_spdf_df) <- c("Value", "x", "y")
  # vie_abnormal <- ggplot() + geom_tile(data = transformed_plastics_spdf_df, aes(x=x,y=y, fill=Value))
  # vie_abnormal
  # # Something's not quite right here, will come back to this later
  # ggarrange(vie_normal, vie_abnormal, nrow = 1, ncol = 2)
  # 
  # sub_raster <- projectRaster(sub_raster, crs = crs(plastics))
  # vie_normal <- ggplot() + geom_raster(data = plastics_spdf_df, aes(x=x,y=y, fill=Value))
  # vie_normal
  # 
  # # Trying out a ggplot
  # # First converting rasters into dataframes- isn't working out
  # sub_raster_df <- as.data.frame(sub_raster, xy = T)
  # transformed_plastics_df <- as.data.frame(transformed_plastics, xy = T)
  # g_trial <- ggplot() + geom_raster(data = transformed_plastics_df, aes(x=x,y=y)) +
  #   scale_fill_gradientn(name = "Floating plastic debris", colors = terrain.colors(10)) +
  #   coord_quickmap()
  #   # geom_raster(data = sub_raster_df, aes(x=x, y=y))
  # g_trial
  # 
  # # Overlaying----
  # sub_raster <- terra::resample(sub_raster, plastics)
  # View(as.data.frame(sub_raster))
  # ggplot() + geom_tile(data = as.data.frame(sub_raster), aes(x=x,y=y))
  # 
  # sub_raster <- terra::project(sub_raster, crs(plastics), res = res(plastics))
  # 
  # 
  # 
  # 
  # # Just changing the crs doesn't reproject the data in itself; we'd have to use the projectRaster function for reprojecting the actual data
  # # plastics_2@crs <- sub_raster@crs
  # 
  # # sub_extent <- as(extent(as.vector(sub_spdf@bbox)), 'SpatialPolygons')
  # # crs(sub_extent) <- sub_crs
  # # plastics_trimmed <- crop(plastics, sub_extent)
  # # str(plastics)
  # # plastics_sp <- as(plastics, "SpatialPixels")
  # # if h was set to 200 km, I get very vague graphs; maybe the null grid settings are wrong; I need to cut the plastics raster precisely according to the extents or bounding boxes set by the kernels to match them up. 
  # 
  
  png(paste0('renditions_output/loop_outputs/kernel_image_',i,'.png'), width = 1379, height = 750)
  plot(unlist(image(kernel)))
  ud <- getverticeshr(kernel, percent = 95)
  ud_sf <- st_as_sf(ud)
  png(paste0('renditions_output/loop_outputs/UD_image_',i,'.png'), width = 1379, height = 750)
  plot(st_geometry(ud_sf[1,]))
  dev.off()
  ud@data$id <- as.factor(ud@data$id)
  png(paste0('renditions_output/loop_outputs/coloured_UD_image_',i,'.png'), width = 1379, height = 750)
  plot(ud, col = ud@data$id)
  dev.off()
  
  # png(paste0('renditions_output/loop_outputs/comparison_between_MCP_and_KDE_',i,'.png'), width = 1379, height = 750)
  # plot(ud, col = "red")
  # plot(mcp(sub_SpatialPoints, percent = 100), add = T)
  # plot(sub_SpatialPoints, cex = 0.75, pch = 1, add = T)
  # dev.off()
}
