# Rescaling plastics raster and revised KDE script

# Install essential packages

library(raster)
library(cowplot)
library(stringr)
library(RColorBrewer)
library(viridisLite)
library(viridis)
library(spData)

# Defining directories
wd <- setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels")
dir_1by1 <- ("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/renditions_output/revised_script_12_1/")

# Loading in relevant_new_data_2----
setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/input_data/fulmars_project_data/")
new_data_1 <- readRDS("test_2colonies.rds")
new_data_2 <- readRDS("test_2colonies_individ_info.rds")
indiv_merged_df <- merge(new_data_1, new_data_2, by = "individ_id")
relevant_new_data_1 <- dplyr::select(indiv_merged_df, individ_id, timestamp, lon, lat, loc_type, colony)
df <- st_as_sf(relevant_new_data_1, coords = c('lon','lat'), crs = 4326)
relevant_new_data_2 <- relevant_new_data_1 %>% filter(!grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp))

# Loading in world 
land <- as(world, "Spatial")

# Loading plastics raster
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
p_sum1 <- plastics2/sum(raster::getValues(plastics2)) # why, to make the values more generalizable? What if I don't do this? 
p_sum1[is.na(plastics)] <- NA

# # Restoring base plotting settings
# old.par <- par(mar = c(0, 0, 0, 0))
# par(old.par)

res(p_sum1)

par(mfrow = c(1,1))
plot(p_sum1)

RES <- res(plastics) # the resolution of the raster (in degrees)

R <- 6371007.2 # the Earth's authalic radius (in meters)- what is authalic???????????????????????????????????????????????????????????????
lat <- raster::yFromRow(plastics, 1:nrow(plastics)) # latitude of the centroid of each cell (in degrees, need to be converted in radians)
area <- (sin(pi/180*(lat + RES[2]/2)) - sin(pi/180*(lat - RES[2]/2))) * (RES[1] * pi/180) * R^2
r_area <- raster::setValues(plastics, rep(area, each = ncol(plastics))) # gives the area of each grid cell in meters 
plot(r_area, col = colsviri)

dir_demClasses <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/input_data/fulmars_project_data/colony_raster_files"
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
  # b <- median(raster::getValues(a)) 
  a_proj <- raster::projectRaster(a, plastics, method = "bilinear")
  a_proj2 <- a_proj * r_area / 100000000 # rescaling the values; dividing by 10^8 because the original raster size was 10 km^2
  a_proj2[is.na(a_proj2)] <- 0 
  
  sum(raster::getValues(a_proj2))
  a_proj2 <- a_proj2/sum(raster::getValues(a_proj2))
    
  plot(a_proj2)
  # c <- median(raster::getValues(a_proj2)) 
  
  # #find number of NAs for plastics where birds are
  # na_a_proj2 <- a_proj2
  # na_a_proj2[na_a_proj2 > 0] <- 1
  # na_a_proj2[is.na(na_a_proj2)] <- 0
  # nas$vals[i] <- sum(na_a_proj2@data@values)
  # 
  # na_p_sum1 <- p_sum1
  # na_p_sum1[is.na(na_p_sum1)] <- 1
  # na_p_sum1[na_p_sum1 != 1] <- 0
  # na_over <- na_a_proj2 * na_p_sum1
  # nas$nas[i] <- sum(na_over@data@values) 
  # 
  # a_proj2[is.na(plastics)] <- NA 
  # Aforementioned part commented out for the time being
  
  ## exporting results
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
  


 
 

