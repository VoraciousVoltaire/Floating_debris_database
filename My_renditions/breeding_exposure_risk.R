# My rendition of the plastic dataset code----

#Import plastics data, average between three models
#Convert into an Atlantic-centred plastic density raster
#Win Cowger & Beth Clark 2021
rm(list=ls()) 

# Load packages and data ####
library(sp)
library(raster)
# library(rgdal) 
library(dplyr)
library(RColorBrewer)
library(sf)

sessionInfo()
#R version 4.1.2 (2021-11-01)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 19045)

#Matrix products: default

#locale:
#[1] LC_COLLATE=English_United Kingdom.1252 
#[2] LC_CTYPE=English_United Kingdom.1252   
#[3] LC_MONETARY=English_United Kingdom.1252
#[4] LC_NUMERIC=C                           
#[5] LC_TIME=English_United Kingdom.1252    

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods  
#[7] base     

#other attached packages:
# [1] RColorBrewer_1.1-2 dplyr_1.0.8        rgdal_1.4-8       
#[4] raster_3.1-5       sp_1.5-0          

#loaded via a namespace (and not attached):
#[1] Rcpp_1.0.8       rstudioapi_0.13  magrittr_2.0.2   tidyselect_1.1.2
#[5] lattice_0.20-45  R6_2.5.1         rlang_1.0.6      fansi_1.0.2     
#[9] tools_4.1.2      grid_4.1.2       utf8_1.2.2       cli_3.3.0       
#[13] DBI_1.1.2        ellipsis_0.3.2   assertthat_0.2.1 tibble_3.1.6    
#[17] lifecycle_1.0.3  crayon_1.5.0     purrr_0.3.4      vctrs_0.3.8     
#[21] codetools_0.2-18 glue_1.6.2       compiler_4.1.2   pillar_1.7.0    
#[25] generics_0.1.2   pkgconfig_2.0.3 

## specify/create directories
dir.create("outputs/")

#Data to read in ----
Lebreton <- as.matrix(read.csv("input_data/plastics_data/lebretonmodel_abundance.csv", header = F))
Maximenko <- as.matrix(read.csv("input_data/plastics_data/maximenkomodel_abundance.csv", header = F))
VanSeb <- as.matrix(read.csv("input_data/plastics_data/vansebillemodel_abundance.csv", header = F))

#Data Cleanup ----
df <- data.frame(van = as.vector(VanSeb), max = as.vector(Maximenko), leb = as.vector(Lebreton))

#Geomean

# My addition 1---- removed 10 to the power in the object Average
Average <- rowMeans(mutate_all(df, function(x) log10(x+1)) ,na.rm = T)

dim(Average) <- c(181, 361)

Ave <- raster(Average, xmn = 1, xmx= 361, ymn=-90, ymx=90, 
              crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

AveRaster_1 <- raster :: shift(rotate(raster(Average, 
                                             xmn = 0, xmx= 360, ymn=-90, ymx=90, 
                                             crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")), dx=0.5)

AveRaster_2 <- raster :: rotate(raster(Average, 
                                       xmn = 0, xmx= 360, ymn=-90, ymx=90, 
                                       crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# My addition 2---- removed log
plot(Ave)
plot(AveRaster_1) # isn't required is my point
plot(AveRaster_2)

#After rotating to atlantic-centred, there is
#a column in the centre with incorrect values
#corrected below

plastics <- AveRaster_1

r178 <- plastics[cellFromCol(plastics,178)]
cols <- as.data.frame(r178)
cols$r179 <- plastics[cellFromCol(plastics,179)]
cols$r180 <- plastics[cellFromCol(plastics,180)]

cols$r182 <- plastics[cellFromCol(plastics,182)]
cols$r183 <- plastics[cellFromCol(plastics,183)]
cols$r184 <- plastics[cellFromCol(plastics,184)]

cols$mean <- rowMeans(cols,na.rm = T)

cols$mean <- ifelse(cols$mean == "NaN",NA,cols$mean)
cols$mean <- ifelse(is.na(cols$r180) & is.na(cols$r182),NA,cols$mean)

plastics[cellFromCol(plastics,181)] <- cols$mean

plot(plastics)

#save the raster
raster_name <- "outputs/00_PlasticsRaster.tif"
writeRaster(plastics, filename = raster_name,
            format="GTiff", overwrite=TRUE)

#Plot difference in coverage between the three models ####

#Read in land file for visualisation:
#Natural Earth land 1:10m polygons version 5.1.1 
#downloaded from www.naturalearthdata.com/

# My addition 3---- rgeoda isn't available in my current version of R version 4.3.1
library(spData)
world_sp <- as(world,"Spatial")

# My addition 4----
# # Isn't required?
# land <- rgdal::readOGR(dsn = "input_data/baselayer", layer = "ne_10m_land") 
# 
# land <- readOGR(dsn = "input_data/baselayer", layer = "ne_10m_land") 


VanSeb_01 <- ifelse(VanSeb>0,1,0)
VanSeb_01 <- ifelse(is.na(VanSeb_01),0,VanSeb_01)

Lebreton_01 <- ifelse(Lebreton>0,1,0)
Lebreton_01 <- ifelse(is.na(Lebreton_01),0,Lebreton_01)

Maximenko_01 <- ifelse(Maximenko>0,1,0)
Maximenko_01 <- ifelse(is.na(Maximenko_01),0,Maximenko_01)

sum01 <- VanSeb_01+Lebreton_01+Maximenko_01
sum01_r <- shift(rotate(raster(sum01, 
                               xmn = 0, xmx= 360, ymn=-90, ymx=90, 
                               crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")), 
                 dx=0.5)

sum01_r <- raster :: shift(rotate(raster(sum01, 
                                         xmn = 0, xmx= 360, ymn=-90, ymx=90, 
                                         crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")), 
                           by=0.5)

yelblus <- c(brewer.pal(n = 5, name = "YlGnBu"),"#00172e")

cols_nmods <- c("#F2F3F400",yelblus[3:5])


plot(sum01_r,col=cols_nmods,breaks = c(-1:3))

# My addition 5----
plot(world_sp, col="grey75", add=T)

png("outputs/00_plastics_model_coverage.png", 
    width=1379,height=750)
plot(sum01_r, col=cols_nmods, breaks = c(-1:3))
plot(land, col="grey75", add=T)
dev.off()


# To add: latitudes and longitudes for each colony; set up a radius - 
# first literature survey to find out colony-specific breeding radii 
# set up specific rasters- find out overlapping regions between the two
# find some methods of aggregations - geometric mean again sounds good 

library(spData)
library(sf)

# Loading in world shape file

world_sp <- as(world,"Spatial")
world_sf <- st_as_sf(world_sp)
plot(world_sp, col = "grey")
  

# Efficient script----

install.packages("googlesheets4")
library(googlesheets4)
refined_datasheet <- read_sheet("https://docs.google.com/spreadsheets/d/1bVSxqMkHXXxFcxjekMXMv5nfoKfJ5yn6y5UztVNUgJY/edit#gid=0")
medlat = median(refined_datasheet$Latitude)
medlon = median(refined_datasheet$Longitude)
proj.laea = paste("+proj=laea +lat_0=",round(medlat), " +lon_0=",round(medlon)," +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km ", sep="")
library(sf)
spatial_obj <- st_as_sf(refined_datasheet, coords = c("Longitude", "Latitude"), crs = 4326)
spatial_laea <- st_transform(spatial_obj, crs = crs(proj.laea))

# Loading in coordinates of colonies
#   faroe_coords_pts <- st_point(c(6.798,61.950))
#   jan_coords_pts <- st_point(c(8.718, 70.921))
#   bear_coords_pts <- st_point(c(-18.956, 74.503))
#   skja_coords_pts <- st_point(c(17.410, 65.990))
#   eyne_coords_pts <- st_point(c(3.115, 59.142))  
#   alke_coords_pts <- st_point(c(-18.459, 79.585))  
#   inis_coords_pts <- st_point(c(10.204, 54.128))  
#   litt_coords_pts <- st_point(c(6.585, 52.136))  
#   isle_coords_pts <- st_point(c(6.559, 57.059))  
#   jars_coords_pts <- st_point(c(5.174, 59.150))  
#   lang_coords_pts <- st_point(c(14.650, 66.350))  
# 
# # Setting up geometry for each colony
#   
#   faroe_geometry <- st_sfc(faroe_coords_pts, crs = "EPSG:4326")
#   jan_geometry <- st_sfc(jan_coords_pts, crs = "EPSG:4326")
#   bear_geometry <- st_sfc(bear_coords_pts, crs = "EPSG:4326")
#   skja_geometry <- st_sfc(skja_coords_pts, crs = "EPSG:4326")
#   eyne_geometry <- st_sfc(eyne_coords_pts, crs = "EPSG:4326")
#   alke_geometry <- st_sfc(alke_coords_pts, crs = "EPSG:4326")
#   inis_geometry <- st_sfc(inis_coords_pts, crs = "EPSG:4326")
#   litt_geometry <- st_sfc(litt_coords_pts, crs = "EPSG:4326")
#   isle_geometry <- st_sfc(isle_coords_pts, crs = "EPSG:4326")
#   jars_geometry <- st_sfc(jars_coords_pts, crs = "EPSG:4326")
#   lang_geometry <- st_sfc(lang_coords_pts, crs = "EPSG:4326")
#   
# # Setting up non-geographical aspect for each colony
#   
#   faroe_nongeographical_attributes <- data.frame(name = "Faroe Islands")
#   jan_nongeographical_attributes <- data.frame(name = "Jan Mayen")  
#   bear_nongeographical_attributes <- data.frame(name = "Bear Island")
#   skja_nongeographical_attributes <- data.frame(name = "Skjalfandi")                                                  
#   eyne_nongeographical_attributes <- data.frame(name = "Eynhallow")      
#   alke_nongeographical_attributes <- data.frame(name = "Alkefjellet")
#   inis_nongeographical_attributes <- data.frame(name = "Inishkea")  
#   litt_nongeographical_attributes <- data.frame(name = "Little Saltee")  
#   isle_nongeographical_attributes <- data.frame(name = "Isle of Canna")  
#   jars_nongeographical_attributes <- data.frame(name = "Jarsteinen")  
#   lang_nongeographical_attributes <- data.frame(name = "Langanes")  
# 
# # Creating an sf object
#   
#   faroe_sf <- st_sf(faroe_nongeographical_attributes, geometry = faroe_geometry)
#   jan_sf <- st_sf(jan_nongeographical_attributes, geometry = jan_geometry)
#   bear_sf <- st_sf(bear_nongeographical_attributes, geometry = bear_geometry)
#   skja_sf <- st_sf(skja_nongeographical_attributes, geometry = skja_geometry)
#   eyne_sf <- st_sf(eyne_nongeographical_attributes, geometry = eyne_geometry)
#   alke_sf <- st_sf(alke_nongeographical_attributes, geometry = alke_geometry)
#   inis_sf <- st_sf(inis_nongeographical_attributes, geometry = inis_geometry)
#   litt_sf <- st_sf(litt_nongeographical_attributes, geometry = litt_geometry)
#   isle_sf <- st_sf(isle_nongeographical_attributes, geometry = isle_geometry)
#   jars_sf <- st_sf(jars_nongeographical_attributes, geometry = jars_geometry)
#   lang_sf <- st_sf(lang_nongeographical_attributes, geometry = lang_geometry)

# Overlaying colonies atop plastic dataset

  # plot(st_geometry(faroe_sf), cex = 0.4, col = "blue", pch = 16, add = T)
  # plot(st_geometry(jan_sf), cex = 0.4, col = "blue", pch = 16, add = T)
  # plot(st_geometry(bear_sf), cex = 0.4, col = "blue", pch = 16, add = T)
  # plot(st_geometry(skja_sf), cex = 0.4, col = "blue", pch = 16, add = T)
  # plot(st_geometry(eyne_sf), cex = 0.4, col = "blue", pch = 16, add = T)
  # plot(st_geometry(alke_sf), cex = 0.4, col = "blue", pch = 16, add = T)
  # plot(st_geometry(inis_sf), cex = 0.4, col = "blue", pch = 16, add = T)
  # plot(st_geometry(litt_sf), cex = 0.4, col = "blue", pch = 16, add = T)
  # plot(st_geometry(isle_sf), cex = 0.4, col = "blue", pch = 16, add = T)
  # plot(st_geometry(jars_sf), cex = 0.4, col = "blue", pch = 16, add = T)
  # overlap_plot <- plot(st_geometry(lang_sf), cex = 0.4, col = "blue", pch = 16, add = T)
  
  raster_name_2 <- "outputs/overlap_plot.tif"
  
# Setting up a buffer of 250 km for each colony
  
  colony_buff <- st_buffer(spatial_laea, dist = 250)
  colony_buff_4326 <- st_transform(colony_buff, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  plot(plastics)
  plot(colony_buff_4326, add = T)
  View(colony_buff_4326)
  library(dplyr)
  colony_buff_refined <- colony_buff_4326 |> select(c('Colony', 'geometry')) 
  
  # Set up a for loop for this
  
  # faroe_crop <- mean(values(crop(plastics, colony_buff_refined[1,2])), na.rm = T)
  # jan_crop <- mean(values(crop(plastics, colony_buff_refined[2,2])), na.rm = T)
  # bear_crop <- mean(values(crop(plastics, colony_buff_refined[3,2])), na.rm = T)
  # skja_crop <- mean(values(crop(plastics, colony_buff_refined[4,2])), na.rm = T)
  # eyne_crop <- mean(values(crop(plastics, colony_buff_refined[5,2])), na.rm = T)
  # alke_crop <- mean(values(crop(plastics, colony_buff_refined[6,2])), na.rm = T)
  # inis_crop <- mean(values(crop(plastics, colony_buff_refined[7,2])), na.rm = T)
  # litt_crop <- mean(values(crop(plastics, colony_buff_refined[8,2])), na.rm = T)
  # isle_crop <- mean(values(crop(plastics, colony_buff_refined[9,2])), na.rm = T)
  # jars_crop <- mean(values(crop(plastics, colony_buff_refined[10,2])), na.rm = T)
  # lang_crop <- mean(values(crop(plastics, colony_buff_refined[11,2])), na.rm = T)

 result <- vector("list", 11)
 for(i in 1:11){result[[i]] <- mean(values
                                    (crop
                                      (plastics, colony_buff_refined[i,2]))
                                      , na.rm = T)}
 means_vector <- as.vector(unlist(result))
 analysis_df <- cbind(colony_buff_refined, Plastic_debris_mean = means_vector)

  # faroe_buff <- st_buffer(faroe_sf, dist = 250000)
  # class(faroe_buff)
  # # Do an st_transform here- as in first set median longitude and latitutde for all points of the colony
  # # then make a crs from the medlong and medlat, then reproject sf obj and their buffers in this new crs (Lambert azimuthal equal-area code)
  # # then st_transform them back to WGS84 (km to degrees conversion again)
  # jan_buff <- st_buffer(jan_sf, dist = 250000)
  # bear_buff <- st_buffer(bear_sf, dist = 250000)
  # skja_buff <- st_buffer(skja_sf, dist = 250000)
  # eyne_buff <- st_buffer(eyne_sf, dist = 250000)
  # alke_buff <- st_buffer(alke_sf, dist = 250000)
  # inis_buff <- st_buffer(inis_sf, dist = 250000)
  # litt_buff <- st_buffer(litt_sf, dist = 250000)
  # isle_buff <- st_buffer(isle_sf, dist = 250000)
  # jars_buff <- st_buffer(jars_sf, dist = 250000)
  # lang_buff <- st_buffer(lang_sf, dist = 250000)
  
# Replacing overlap_plot with buff_overlap_plot

  # plot(st_geometry(faroe_buff), col = "blue", add = T)
  # plot(st_geometry(jan_buff), col = "blue", add = T)
  # plot(st_geometry(bear_buff), col = "blue", add = T)
  # plot(st_geometry(skja_buff), col = "blue", add = T)
  # plot(st_geometry(eyne_buff), col = "blue", add = T)
  # plot(st_geometry(alke_buff), col = "blue", add = T)
  # plot(st_geometry(inis_buff), col = "blue", add = T)
  # plot(st_geometry(litt_buff), col = "blue", add = T)
  # plot(st_geometry(isle_buff), col = "blue", add = T)
  # plot(st_geometry(jars_buff), col = "blue", add = T)
  # buff_overlap_plot <- plot(st_geometry(lang_buff), col = "blue", add = T)
  
# Cropping colony-specific rasters from plastics dataset

#   mean_faroe_cropped <- mean(values(crop(plastics, buff[1,])), na.rm = T)
#   mean_jan_cropped<- mean(values(crop(plastics, buff[2,])), na.rm = T)
#   mean_bear_cropped <- mean(values(crop(plastics, buff[3,])), na.rm = T)
#   mean_skja_cropped <- mean(values(crop(plastics, buff[4,])), na.rm = T)
#   mean_eyne_cropped <- mean(values(crop(plastics, buff[5,])), na.rm = T)
#   mean_alke_cropped <- mean(values(crop(plastics, buff[6,])), na.rm = T)
#   mean_inis_cropped <- mean(values(crop(plastics, buff[7,])), na.rm = T) 
#   mean_litt_cropped <- mean(values(crop(plastics, buff[8,])), na.rm = T)
#   mean_isle_cropped <- mean(values(crop(plastics, buff[9,])), na.rm = T)
#   mean_jars_cropped <- mean(values(crop(plastics, buff[10,])), na.rm = T) 
#   mean_lang_cropped <- mean(values(crop(plastics, buff[11,])), na.rm = T)
#   
# faroe_cropped <- crop(plastics, faroe_buff)
# View(faroe_buff)

#   jan_cropped <- crop(plastics, jan_buff)
#   bear_cropped <- crop(plastics, bear_buff)
#   skja_cropped <- crop(plastics, skja_buff)
#   eyne_cropped <- crop(plastics, eyne_buff)
#   alke_cropped <- crop(plastics, alke_buff)
#   inis_cropped <- crop(plastics, inis_buff)  
#   litt_cropped <- crop(plastics, litt_buff)  
#   isle_cropped <- crop(plastics, isle_buff)  
#   jars_cropped <- crop(plastics, jars_buff)  
#   lang_cropped <- crop(plastics, lang_buff)  
#   some_vector <- c(faroe_cropped, jan_cropped, bear_cropped, skja_cropped,
#                    eyne_cropped, alke_cropped, inis_cropped, litt_cropped, 
#                    isle_cropped, jars_cropped, lang_cropped)
# summary(some_vector)

# Creating a values' vector
  # value_1 <- values(faroe_cropped) 
  # value_2 <- values(jan_cropped) 
  # value_3 <- values(bear_cropped) 
  # value_4 <- values(skja_cropped) 
  # value_5 <- values(eyne_cropped) 
  # value_6 <- values(alke_cropped) 
  # value_7 <- values(inis_cropped) 
  # value_8 <- values(litt_cropped) 
  # value_9 <- values(isle_cropped) 
  # value_10 <- values(jars_cropped) 
  # value_11 <- values(lang_cropped) 
  # values_list <- list(c(value_1, value_2, value_3, value_4, value_5, value_6,
  #                    value_7, value_8, value_9, value_10, value_11))
  
# Calculating means of values inside each colony-specific raster
  
  # mean_faroe_cropped <- mean(values(faroe_cropped), na.rm = T)
  # mean_jan_cropped <- mean(values(jan_cropped), na.rm = T)
  # mean_bear_cropped <- mean(values(bear_cropped), na.rm = T)
  # mean_skja_cropped <- mean(values(skja_cropped), na.rm = T)
  # mean_eyne_cropped <- mean(values(eyne_cropped), na.rm = T)
  # mean_alke_cropped <- mean(values(alke_cropped), na.rm = T)
  # mean_inis_cropped <- mean(values(inis_cropped), na.rm = T)
  # mean_litt_cropped <- mean(values(litt_cropped), na.rm = T)
  # mean_isle_cropped <- mean(values(isle_cropped), na.rm = T)
  # mean_jars_cropped <- mean(values(jars_cropped), na.rm = T)
  # mean_lang_cropped <- mean(values(lang_cropped), na.rm = T)

  # analysis_df <- data.frame(Colonies = as.factor(c("Faroe islands", "Jan Mayen", "Bear island", 
  #                                        "Skjalfandi", "Eynehallow", "Alkefjellet",
  #                                        "Inishkea", "Little Saltee", "Isle of Canna",
  #                                        "Jarsteinen", "Langanes")), 
  #                           Plastic_debris_mean = c(mean_faroe_cropped,  mean_jan_cropped,
  #                                                   mean_bear_cropped, mean_skja_cropped,
  #                                                   mean_eyne_cropped, mean_alke_cropped,
  #                                                   mean_inis_cropped, mean_litt_cropped,
  #                                                   mean_isle_cropped, mean_jars_cropped,
  #                                                   mean_lang_cropped))
 
  # Ordering colony-wise exposure in an intact order
  analysis_df$Colony <- factor(analysis_df$Colony, levels = analysis_df$Colony)
  str(analysis_df)

# Notes ----
# Less coverage for Bear islands, Alkefjellet

# Plotting a basic ggplot for analysis_df

library(ggplot2)
analysis_df_bar_plot <- ggplot(data = analysis_df, aes(x = Colony, y = Plastic_debris_mean)) +
  geom_col(colour = "black", fill = "skyblue") +
  scale_y_continuous(name = "Plastic debris mean") +
  coord_cartesian(ylim = c(1.1,5)) +
  theme_classic() +
  theme(axis.line = element_line(colour = "grey"), axis.text.x = element_text(angle = 90)) 
analysis_df_bar_plot

# Defining a standard deviation column in analysis_df for error bars in box-plot

# sd1 <- sd(values(faroe_cropped), na.rm = T)
# sd2 <- sd(values(jan_cropped), na.rm = T)
# sd3 <- sd(values(bear_cropped), na.rm = T)
# sd4 <- sd(values(skja_cropped), na.rm = T)
# sd5 <- sd(values(eyne_cropped), na.rm = T)
# sd6 <- sd(values(alke_cropped), na.rm = T)
# sd7 <- sd(values(inis_cropped), na.rm = T)
# sd8 <- sd(values(litt_cropped), na.rm = T)
# sd9 <- sd(values(isle_cropped), na.rm = T)
# sd10 <- sd(values(jars_cropped), na.rm = T)
# sd11 <- sd(values(lang_cropped), na.rm = T)
# analysis_df_2 <- cbind(analysis_df, Standard_devation = c(sd1, sd2, sd3, sd4, sd5, sd6, sd7, 
#                                                           sd8, sd9, sd10, sd11))

sd_result <- vector("list", 11)
for(i in 1:11){sd_result[[i]] <- sd(values
                                   (
                                     crop
                                     (plastics, colony_buff_refined[i,2])
                                     )
                                   , na.rm = T)}
sd_vector <- as.vector(unlist(sd_result))
analysis_df_2 <- cbind(analysis_df, Standard_deviation = sd_vector)
transmute_seom <- analysis_df_2 |> transmute(Standard_error_of_mean = Standard_deviation/sqrt(Sample_size))
analysis_df_2 <- cbind(analysis_df_2, transmute_seom)
View(analysis_df_2)

# Trying out a 'box-plot' with standard deviation
  analysis_df_box_plot <- ggplot(data = analysis_df_2, aes(x = Colony, y = Plastic_debris_mean,
                               ymin = Plastic_debris_mean - Standard_deviation,
                               ymax = Plastic_debris_mean + Standard_deviation)) +
  scale_y_continuous(name = "Plastic debris mean") +
  geom_boxplot() +
  theme_classic() +
  theme(axis.line = element_line(colour = "grey"), axis.text.x = element_text(angle = 90)) +
  geom_errorbar()
analysis_df_box_plot

# Trying out a plot with seom
analysis_df_seom <- analysis_df_box_plot <- ggplot(data = analysis_df_2, aes(x = Colony, y = Plastic_debris_mean,
                                                                             ymin = Plastic_debris_mean - Standard_error_of_mean,
                                                                             ymax = Plastic_debris_mean + Standard_error_of_mean)) +
  scale_y_continuous(name = "Plastic debris mean") +
  geom_boxplot() +
  theme_classic() +
  theme(axis.line = element_line(colour = "grey"), axis.text.x = element_text(angle = 90)) +
  geom_errorbar()
analysis_df_seom
 
# Analysis----

# Checking sample size 

# faroe_sum <- sum(!is.na(values(faroe_cropped))) # 30
# jan_sum <- sum(!is.na(values(jan_cropped))) # 56
# bear_sum <- sum(!is.na(values(bear_cropped))) # 24
# skja_sum <- sum(!is.na(values(skja_cropped))) # 18
# eyne_sum <- sum(!is.na(values(eyne_cropped))) # 31
# alke_sum <- sum(!is.na(values(alke_cropped))) # 5
# inis_sum <- sum(!is.na(values(inis_cropped))) # 26
# litt_sum <- sum(!is.na(values(litt_cropped))) # 20
# isle_sum <- sum(!is.na(values(isle_cropped))) # 26
# jars_sum <- sum(!is.na(values(jars_cropped))) # 24
# lang_sum <- sum(!is.na(values(lang_cropped))) # 40

ss_result <- vector("list", 11)
for(i in 1:11){
  ss_result[[i]] <- sum(!is.na
                   (values(crop(plastics, colony_buff_refined[i,2])))
                     )
}
Sample_size <- as.vector(unlist(ss_result))
analysis_df_2 <- cbind(analysis_df, Sample_size)
View(analysis_df_2)

# # Adding sample size to analysis_df
# analysis_df <- cbind(analysis_df, Sample_size = c(faroe_sum, jan_sum, bear_sum, skja_sum, eyne_sum,
#                                     alke_sum, inis_sum, litt_sum, isle_sum, jars_sum,
#                                     lang_sum))
# View(analysis_df)

# Checking the normality assumption for each colony----

# shapiro.test(values(faroe_cropped))
# shapiro.test(values(jan_cropped))
# shapiro.test(values(bear_cropped))
# shapiro.test(values(skja_cropped))
# shapiro.test(values(eyne_cropped)) # Normally distributed
# shapiro.test(values(alke_cropped)) # Normally distributed
# shapiro.test(values(inis_cropped)) # Normally distributed
# shapiro.test(values(litt_cropped)) # Normally distributed
# shapiro.test(values(isle_cropped)) # Normally distributed- but variables don't have to be normally distributed 
#                                    # for running parametric tests, their residuals have to be
# shapiro.test(values(jars_cropped))
# shapiro.test(values(lang_cropped))

library(ggpubr)
dplot <- vector('list', 11)
 for(i in 1:11){ 
   dplot[[i]] <- local({
     i <- i
     print(ggdensity(na.omit(values(crop(plastics, colony_buff_refined[i,2])))), 
                     main = "Density plot",
                     xlab = "Plastic debris value")})
 }
library(gridExtra)
for(j in 1:11){
  grid.arrange(grobs = dplot, ncol = 5, nrow = 3)
}
         
# plot_2 <- ggdensity(values(crop(plastics, colony_buff_refined[2,2])), 
#                     main = "Density plot",
#                     xlab = "Plastic debris value")
# plot_3 <- ggdensity(values(crop(plastics, colony_buff_refined[3,2])), 
#                     main = "Density plot",
#                     xlab = "Plastic debris value")
# plot_4 <- ggdensity(values(crop(plastics, colony_buff_refined[4,2])), 
#                     main = "Density plot",
#                     xlab = "Plastic debris value")
# plot_5 <- ggdensity(values(crop(plastics, colony_buff_refined[5,2])), 
#                     main = "Density plot",
#                     xlab = "Plastic debris value")
# plot_6 <- ggdensity(values(crop(plastics, colony_buff_refined[6,2])), 
#                     main = "Density plot",
#                     xlab = "Plastic debris value")
# plot_7 <- ggdensity(values(crop(plastics, colony_buff_refined[7,2])), 
#                     main = "Density plot",
#                     xlab = "Plastic debris value")
# plot_8 <- ggdensity(values(crop(plastics, colony_buff_refined[8,2])), 
#                     main = "Density plot",
#                     xlab = "Plastic debris value")
# plot_9 <- ggdensity(na.omit(values(crop(plastics, colony_buff_refined[9,2]))), 
#                     main = "Density plot",
#                     xlab = "Plastic debris value")
# plot_10 <- ggdensity(values(crop(plastics, colony_buff_refined[10,2])), 
#                     main = "Density plot",
#                     xlab = "Plastic debris value")
# plot_11 <- ggdensity(values(crop(plastics, colony_buff_refined[11,2])), 
#                     main = "Density plot",
#                     xlab = "Plastic debris value")

shapiro_result <- vector("list", 11)
for(i in 1:11){
  shapiro_result[[i]] <- shapiro.test(values(crop(plastics, colony_buff_refined[i,2])))
}
for(i in 1:11){print(shapiro_result[[i]])}

# Results:
# 1: Non-normal
# 2: Non-normal
# 3: Non-normal
# 4: Non-normal
# 5: Normal
# 6: Normal
# 7: Normal
# 8: Normal
# 9: Normal
# 10: Normal
# 11: Non-normal

# QQ plots ----
library(ggpubr)
library(ggplot2)
qplot <- vector('list', 11)
for(i in 1:11){ 
  qplot[[i]] <- local
  (
    {
       i <- i
       print(ggqqplot(na.omit(values(crop(plastics, colony_buff_refined[i,2])))), 
       main = "Density plot")
    }
  )
}
library(gridExtra)
for(j in 1:11){
  grid.arrange(grobs = qplot, ncol = 5, nrow = 3)
}

# Creating a new analysis data frame with short code
library(terra)
values_list <- vector('list', 269)
for(i in 1:11){
  values_list[[i]] <- na.omit(values(crop(plastics, colony_buff_refined[i,2])))
  }
values_vector <- unlist(values_list)
new_analysis_df <- data.frame(Colony = rep(analysis_df_2$Colony, analysis_df_2$Sample_size),
                              Plastic_debris_value = values_vector)
View(new_analysis_df)

# Performing the non-parametric Kruskal Wallis test
# First create a new dataframe that contains two columns: plastic abundance value
# and colony name (like group)

# faroe_value_df <- na.omit(data.frame(Colony_name = rep("Faroe islands", length(values(faroe_cropped))), 
#                                      Values = values(faroe_cropped), Sample_size = 30))
# jan_value_df <- na.omit(data.frame(Colony_name = rep("Jan Mayen", length(values(jan_cropped))), 
#                                      Values = values(jan_cropped), Sample_size = 56))
# bear_value_df <- na.omit(data.frame(Colony_name = rep("Bear island", length(values(bear_cropped))), 
#                                      Values = values(bear_cropped), Sample_size = 24))
# skja_value_df <- na.omit(data.frame(Colony_name = rep("Skjalfandi", length(values(skja_cropped))), 
#                                      Values = values(skja_cropped), Sample_size = 18))
# eyne_value_df <- na.omit(data.frame(Colony_name = rep("Eynehallow", length(values(eyne_cropped))), 
#                                      Values = values(eyne_cropped), Sample_size = 31))
# alke_value_df <- na.omit(data.frame(Colony_name = rep("Alkefjellet", length(values(alke_cropped))), 
#                                      Values = values(alke_cropped), Sample_size = 5))
# inis_value_df <- na.omit(data.frame(Colony_name = rep("Inishkea", length(values(inis_cropped))), 
#                                      Values = values(inis_cropped), Sample_size = 26))
# litt_value_df <- na.omit(data.frame(Colony_name = rep("Little Saltee", length(values(litt_cropped))), 
#                                      Values = values(litt_cropped), Sample_size = 20))
# isle_value_df <- na.omit(data.frame(Colony_name = rep("Isle of Canna", length(values(isle_cropped))), 
#                                      Values = values(isle_cropped), Sample_size = 26))
# jars_value_df <- na.omit(data.frame(Colony_name = rep("Jarsteinen", length(values(jars_cropped))), 
#                                      Values = values(jars_cropped), Sample_size = 24))
# lang_value_df <- na.omit(data.frame(Colony_name = rep("Langanes", length(values(lang_cropped))), 
#                                      Values = values(lang_cropped), Sample_size = 40))
# 
# merged_value_df <- rbind(faroe_value_df, jan_value_df, bear_value_df, skja_value_df,
#                          eyne_value_df, alke_value_df, inis_value_df, litt_value_df,
#                          isle_value_df, jars_value_df, lang_value_df)
# View(merged_value_df)           

kruskal.test(Plastic_debris_value ~ Colony, data = new_analysis_df)
pairwise.wilcox.test(new_analysis_df$Plastic_debris_value, new_analysis_df$Colony,
                     p.adjust.method = "BH")
plot(glm(Plastic_debris_value ~ Colony, data = new_analysis_df))

kruskal.test(Plastic_debris_value ~ Colony, data = trimmed_new_analysis_df)
pairwise.wilcox.test(trimmed_new_analysis_df$Plastic_debris_value, trimmed_new_analysis_df$ Colony,
                     p.adjust.method = "BH")




library(ggplot2)
boxplot_plastic_debris_values <- ggplot(data = new_analysis_df, aes(x = Colony, y = Plastic_debris_value)) +
  geom_point() +
  geom_boxplot() +
  theme_classic() +
  scale_y_continuous(name = "Plastic debris value") +
  theme(axis.line = element_line(colour = "grey"), axis.text.x = element_text(angle = 90)) 
boxplot_plastic_debris_values


# Results show that Jarsteinen-Eynehallow, Jarsteinen-Faroe islands, 
# Little-Saltee-Inishkea, Skjalfandi-Eynehallow, Skjalfandi-Faroe islands,
# Jarsteinen-Jan Mayen, Langanes-Isle of Canna, Skjalfandi-Isle of Canna,
# Skjalfandi-Jan Mayen, Skjalfandi-Jarsteinen and Skjalfandi-Langanes aren't 
# significantly different; i.e 11 out of 55 total pairs

# Now, performing a non-parametric equivalent of ANCOVA on the same dataset 
# controlling for sample size to see if the significances hold true. Also
# don't forget to apply multiple-testing correction

install.packages("npsm")
install.packages("Rfit")
install.packages("fANCOVA")

# Plot to tease out something- would colouring by covariate help see any trends?

library(ggplot2)
covariate_plot <- ggplot(merged_value_df, aes(x = Colony_name, y = Values, colour = Sample_size)) + 
  geom_point() +
  theme_classic() +
  xlab("Colony name") +
  ylab("Plastic debris values") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Covariate trial", subtitle = "nope")
 covariate_plot 

# Assuming what I did was wrong so far (interpreting normality in variables rather than the residuals)
# I'm now trying to plot residual plots from a generalized linear model 
 
trial_model <- lm(Values ~ Colony_name, data = merged_value_df)
res <- resid(trial_model)
plot(fitted(trial_model), res)
abline(0,0)
qqnorm(res)  
qqline(res)  
plot(density(res))

# This shows that parametric tests can be applied here? I need more reading on this. I've learnt
# not to associate non-parametric and parametric with non-normal and normal respectively but I'm
# having a hard time interpreting tests. They say non-parametric tests sometimes have more assumptions than
# parametric tests? I used to think that just checking normality of data (residuals more technically,
# but I used to think that Shapiro-Wilk test takes care of that and directly tells us whether your 
# data is normal or not based on its residuals- have to check exactly what it is that it does), 
# no heteroskedasticity (Levene's test) and independence of independent variables (sounds funny)
# were the only three assumptions for parametric tests and that if your data violates even one of 'em
# it's straight down to non-parametric tests. Turns out, no tis not so easy and minor violations (depends 
# on the category violated) are permissible to some extent. So I'm trying parametric tests now because
# the density plot clearly resembled a bell-shaped curve

# ANCOVA assumptions: 1) Independence of treatment (plastic debris abundance- proportional to
# exposure- in this case, and covariate- sample size in this case; 2) Homogenity of variance 

# Second assumption testing
library(car)
leveneTest(Values~Colony_name, data = merged_value_df)
# Thus even after transforming (plastic debris abundance data is log transformed), ANCOVA's 
# second assumption is violated; hence doing the next part just for practice

ancova_model <- aov(Values ~ Colony_name + Sample_size, data = merged_value_df)
Anova(ancova_model) # didn't mention the type ("III" in the example) because I don't know what it does
# Controlling for sample size, values and colony name still show a statistically significant
# relation, hence different colonies have different plastic exposures never mind their sample
# sizes but ancova can't be performed here technically because its second assumption is violated.
# I don't know what I'll do at this point, need to read some basics for now.

# I don't need an ANCOVA, rather I needed a factorial ANOVA because all my IVs were categorical;
# however, I have moved on to analysis rather than tests. Hence, I shall fit a generalized linear
# model.

# Robustyfying the data
require(MASS)
robust_model = rlm(Plastic_debris_value ~ Colony, data = new_analysis_df)
summary(robust_model)
plot(robust_model)
# Isn't needed actually because the data is already log transformed and robustyfying it along
# with that isn't recommended?

# In order to not use Kruskal-Wallis, I'd have to use a glm that assumes the residuals to follow
# a non-linear distribution: poisson, gamma, binomial, negative binomial distribution- use compare.fits
# or model.comparison 

# Okay, understanding GLMs---- 
# I think I'll spend time on this later because the gamma family isn't really fitting 
# the data very well (from Q-Q residuals), so for now Kruskal-Wallis saves the day
install.packages("DHARMa")
library(DHARMa)
mod1 <- glm(Plastic_debris_value ~ Colony, family=Gamma, data = trimmed_new_analysis_df)
summary(mod1)
plot(mod1)
simulationOutput <- simulateResiduals(fittedModel = mod1)
plot(simulationOutput)

count_mod <- glm(Plastic_debris_value ~ 1, data = new_analysis_df)
summary(count_mod)
library(sf)





