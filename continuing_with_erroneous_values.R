# Combining month-wise exposure risk scores into populations

# Loading essential packages
library(raster)
library(tidyverse)
library(stringr)

# Direction to rasters
dir_1by1 <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/renditions_output/last_resort/fidelity/KDE_loop_2/latest_attempt_tifs_and_pics_both_normalized"
dir_out <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/renditions_output/last_resort/fidelity/month_summed_rasters"
dat <- read.csv("/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/renditions_output/last_resort/fidelity/exposure_scores_by_month.csv")
head(dat)
summary(dat$exposure_score)

# Combine by population

pop_exposure <- dat %>%
  group_by(population) %>%  
  summarise(population_exposure = mean(exposure_score)) %>%
  data.frame() 

pop_exposure$density_sum <- NA
pop_exposure$new_values <- NA

pop_exposure$population <- gsub("\\."," ",pop_exposure$population) # replacing dot with a space
months <- list.files(dir_1by1, pattern = paste0(".*\\.tif$"))


# setwd(paste0(dir_1by1))
# for(p in pop_exposure$population){
#   rasters <- months[grepl(p, months, fixed = T)]
#   raster_sum <-  # creating a raster stack and then doing a mean operation is better
#   raster_name <- paste0(dir_out,"/",pop_exposure$population[i],".tif")
#   raster::writeRaster(rast_sum, filename = raster_name, format="GTiff", overwrite = TRUE)
#   
#   rast_sum[is.na(rast_sum)] <- 0
#   density_sum <- sum(raster::getValues(rast_sum))
#   pop_exposure$density_sum[i] <- density_sum
#   
#   print(pop_exposure$sp_pop[i])
# }


for (i in 1:nrow(dat)){
    months <- list.files(dir_1by1, pattern = paste0(".*\\.tif$"))
}
View(as.data.frame(months))
 
  
# Define p_sum1 here

  # Working on this for loop- takes a lot of time to run, returning mostly empty rasters
for(i in 1:nrow(pop_exposure)){
  
  rasters <- months[grepl(pop_exposure$population[i], months, fixed = T)]
  
  for(j in 1:length(rasters)){
    a <- raster::raster(paste0(dir_1by1,"/", rasters[j]))
    a <- raster::projectRaster(a, p_sum1, method = "bilinear")
    a*r_area/ 100000000
    
    if(j == 1){rast_sum <- a} else{rast_sum <- a + rast_sum}
    
    rast_sum[is.na(rast_sum)] <- 0
    print(i)
    print(max(rast_sum@data@values))
    
    weighted_rast_sum <- rast_sum
    
    # weighted_rast_sum <- rast_sum/6 # Try not taking the mean, we'll deal with Alkefjellet later
    raster_name <- paste0(dir_out,"/",pop_exposure$population[i],".tif")
    raster::writeRaster(weighted_rast_sum, filename = raster_name, format="GTiff", overwrite = TRUE)
    density_sum <- sum(raster::getValues(weighted_rast_sum))
    pop_exposure$density_sum[i] <- density_sum
    new_score_plot <- weighted_rast_sum * p_sum1
    new_score <- new_score_plot
    new_score[is.na(new_score)] <- 0
    
    # I don't believe that the plastic exposure values won't change after doing this; thus multiplying this with the plastics raster
   
    new_values <- round(sum(raster::getValues(new_score))*1000000, 4) # * This method is different how from the last method? In the last method, we normalized both values (didn't really matter then because the denominator was almost 1), multiplies, got an exposure risk score for each month and then added the scores together. In this method, we first summed the rasters and then multiplied them with the plastics raster. Ultimately, the value arrived must be the SAME
    pop_exposure$new_values[i] <- new_values
   
    png(paste0(dir_out, "/", pop_exposure$population[i], ".png"), width = 1399, height = 455)
    par(mfrow=c(1,2))
    plot(weighted_rast_sum, main = paste0(pop_exposure$population[i]," distribution"),col = colsviri,legend = F)
    plot(new_score_plot, main = paste0("Exposure score = ", new_values),
         col = colsinf, legend = F)
    dev.off()
  }
} 

View(pop_exposure)
range(pop_exposure$density_sum)

  
  
  
  
 

  
  print(pop_exposure$sp_pop[i])
  print(i)


range(pop_exposure$density_sum)


