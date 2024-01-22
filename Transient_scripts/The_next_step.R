
# Install necessary packages----

library(raster)
library(terra)
library(RColorBrewer)
library(sf)
library(tidyverse)
library(sp)
library(viridisLite)
library(viridis)
library(cowplot)
library(spData)
library(rnaturalearth)

sessionInfo()

# Loading plastics raster----
setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/")
plastics <- raster("outputs/00_PlasticsRaster.tif")

## Rescale to 1
plastics2 <- plastics
plastics2[is.na(plastics2)] <- 0 # NA turned to 0s
p_sum1    <- plastics2/sum(raster::getValues(plastics2)) # getting quantitative value proportions
p_sum1[is.na(plastics)] <- NA # for any weird divisions

yelblus <- c(brewer.pal(n = 9, name = "YlGnBu"),"#00172e") # from the RColorBrewer package
col_birds <- c(colorRampPalette(yelblus)(1000))
plot(plastics2)

# Loading in world 
land <- as(world, "Spatial")

# Define Robinson projection
proj <- "+proj=robin"

# Loading in relevant_new_data_2----
setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/input_data/fulmars_project_data/")
new_data_1 <- readRDS("test_2colonies.rds")
new_data_2 <- readRDS("test_2colonies_individ_info.rds")
indiv_merged_df <- merge(new_data_1, new_data_2, by = "individ_id")
relevant_new_data_1 <- dplyr::select(indiv_merged_df, individ_id, timestamp, lon, lat, loc_type, colony)
df <- st_as_sf(relevant_new_data_1, coords = c('lon','lat'), crs = 4326)
relevant_new_data_2 <- relevant_new_data_1 %>% filter(!grepl(c('-04-|-05-|-06-|-07-|-08-|-09-') ,timestamp))

# Too advanced for me right now; I'd first have to bring the tif figures that I've got which represent the volume rasters from getvolumeUD
setwd("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/renditions_output/")
for(i in unique(relevant_new_data_2$colony)){
  sub <- raster(paste0("revised_script_12_1/",i,".tif"))
  plot(sub) 
  Sys.sleep(2)
}

par(mfrow=c(1,1))

# Trying a ggplot- ggplot's really messing up these rasters. If I run individual lines of a loop, it's working fine but the computational power is lacking to graph looped outputs it seems.
# Okay, my bad. Just had to take the plastics_dup line inside the loop

j <- 1
npers <- vector("list", length(unique(relevant_new_data_2$colony)))
pers <- vector("list", length(unique(relevant_new_data_2$colony)))
for(i in unique(relevant_new_data_2$colony)){
  # no_of_iterations <- 0
  sub <- raster(paste0("revised_script_12_1/",i,".tif"))
  
  # plastics_dup <- plastics
  # plastics_dup <- projectRaster(plastics_dup, crs = crs(sub))
  sub <- resample(sub, plastics, method = "bilinear")
  sub.df <- as.data.frame(sub, xy = T)
  colnames(sub.df) <- c("x","y","kernelUDvolume")
  
 
  
  plastics_3 = as.data.frame(plastics, xy = T)
  colnames(plastics_3) <- c("x","y","Value")
  sub.df[is.na(sub.df)] <- 0
  print(
    ggplot() + geom_raster(data = plastics_3, aes(x=x, y=y, fill = Value)) + 
      scale_fill_viridis_c() +
      geom_raster(data = sub.df, aes(x=x,y=y, alpha = kernelUDvolume)) +
      geom_sf(data = ne_countries(scale = 50)) + 
          coord_sf(xlim = c(as.vector(sub@extent)[1], 
                            as.vector(sub@extent)[2]), 
                   ylim = c(as.vector(sub@extent)[3], 
                            as.vector(sub@extent)[4]))
  )
  
          

  
  
  # print(ggplot() + geom_raster(data = plastics_dup.df, aes(x=x, y=y, fill = Value)) +
  #         geom_raster(data = sub.df, aes(x=x, y=y, alpha = kernelUDvolume)) +
  #         # scale_fill_gradient(name = "Floating plastics debris concentration", low = "white", high = "red") +
  #         scale_fill_viridis_c() +
  #         geom_sf(data = ne_countries(scale = 50)) +
  #         coord_sf(xlim = c(min(sub.df$x) - 5, max(sub.df$x) + 5), ylim = c(min(sub.df$y) - 5, max(sub.df$y) + 5)) + 
  #         coord_quickmap() +
  #         theme_minimal() +
  #         ggtitle(paste0("Colony: ",i)) + 
  #         theme(axis.title = element_blank()))
  
  
  # no_of_iterations <- no_of_iterations + 1
  # if(no_of_iterations == 1){
  #   break
  # }
  
  # Sys.sleep(2)
  ggsave(paste0('revised_script_12_1/loop_outputs/overlapped_rasters_',i,'.png'))
  
  new_raster <- plastics_dup * sub
  new_raster_df <- as.data.frame(new_raster, xy = T)
  colnames(new_raster_df) <- c("x", "y", "Value")
  print(ggplot() + geom_raster(data = new_raster_df, aes(x=x,y=y, fill = Value)) + 
          geom_sf(data = ne_countries(scale = 50)) + 
          coord_sf(xlim = c(as.vector(new_raster@extent)[1] + 5, 
                            as.vector(new_raster@extent)[2]) - 3, 
                   ylim = c(as.vector(new_raster@extent)[3] + 1, 
                            as.vector(new_raster@extent)[4] - 1)))
  ggsave(paste0('revised_script_12_1/loop_outputs/multiplication_raster_',i,'.png'))
        
  npers[[j]] <- median(new_raster_df$Value, na.rm = T)
  print(npers[[j]])
  # Creating a dataset which stores these multiplicated values to run regression analysis
  pers[[j]] <- as.vector(na.omit(new_raster_df$Value))
  j <- j + 1 
  
  
  
  # print(ggplot() + geom_raster(data = plastics_dup.df, aes(x=x, y=y, fill = Value)) +
  #   geom_raster(data = sub.df, aes(x=x, y=y, alpha = kernelUDvolume)) +
  #   # scale_fill_gradient(name = "Floating plastics debris concentration", low = "white", high = "red") +
  #   scale_fill_viridis_c() +
  #   coord_quickmap() +
  #   theme_minimal() +
  #     ggtitle(paste0("Colony: ",i)) + 
  #   theme(axis.title = element_blank()))
  
  # npers <- vector("list", length(unique(relevant_new_data_2$colony)))
  # for(j in 1:length(unique(relevant_new_data_2$colony))){
  #   if(unique(relevant_new_data_2$colony[j] == i)){
  # npers[[j]] <- median(new_raster_df$Value, na.rm = T)
  #     print(npers[[j]])
  #   }
  # }
}

Pseudo_analysis_df <- data.frame("Colonies" = unique(relevant_new_data_2$colony), "Plastic exposure risk score" = as.vector(unlist(npers)))
View(Pseudo_analysis_df)

# I'm getting erroneous graphs after this
# Next steps: change the crs, ocrrect ggplot extents, merge data frames, figure out hr95





