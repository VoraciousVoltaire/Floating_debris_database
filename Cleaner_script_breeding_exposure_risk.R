# Loading in essential packages ----
library(sp)
library(raster)
library(dplyr)
library(RColorBrewer)
library(sf)
library(spData)
library(googlesheets4)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(terra)

# Setting up the right directory ----
setwd('/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/1_full_analysis_petrels/')

# Data to read in ----
Lebreton <- as.matrix(read.csv("input_data/plastics_data/lebretonmodel_abundance.csv", header = F))
Maximenko <- as.matrix(read.csv("input_data/plastics_data/maximenkomodel_abundance.csv", header = F))
VanSeb <- as.matrix(read.csv("input_data/plastics_data/vansebillemodel_abundance.csv", header = F))

# Data Cleanup ----
df <- data.frame(van = as.vector(VanSeb), max = as.vector(Maximenko), leb = as.vector(Lebreton))

# Geometric mean ----
Average <- rowMeans(mutate_all(df, function(x) log10(x+1)) ,na.rm = T)

dim(Average) <- c(181, 361)

Ave <- raster(Average, xmn = 1, xmx= 361, ymn=-90, ymx=90, 
              crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(Ave)

AveRaster_2 <- raster :: rotate(raster(Average, 
                                       xmn = 0, xmx= 360, ymn=-90, ymx=90, 
                                       crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(AveRaster_2)

# After rotating to atlantic-centred, there is a column in the center with incorrect values
# corrected below ---- 

plastics <- AveRaster_2

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

# Plot difference in coverage between the three models ---- 
world_sp <- as(world,"Spatial")

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
plot(world_sp, col="grey75", add=T)

# Reading in data from refined_datasheet ----

refined_datasheet <- read_sheet("https://docs.google.com/spreadsheets/d/1bVSxqMkHXXxFcxjekMXMv5nfoKfJ5yn6y5UztVNUgJY/edit#gid=0")

medlat = median(refined_datasheet$Latitude)
medlon = median(refined_datasheet$Longitude)

proj.laea = paste("+proj=laea +lat_0=",round(medlat), " +lon_0=",round(medlon)," +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km ", sep="")
spatial_obj <- st_as_sf(refined_datasheet, coords = c("Longitude", "Latitude"), crs = 4326)
spatial_laea <- st_transform(spatial_obj, crs = crs(proj.laea))

# Setting up a buffer of 250 km for each colony ----

colony_buff <- st_buffer(spatial_laea, dist = 250)
colony_buff_4326 <- st_transform(colony_buff, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(plastics)
plot(colony_buff_4326, add = T)
colony_buff_refined <- colony_buff_4326 |> select(c('Colony', 'geometry')) 

# Mean plastic abundance value for each colony ----
result <- vector("list", 11)
for(i in 1:11){result[[i]] <- mean(values
                                   (crop
                                     (plastics, colony_buff_refined[i,2]))
                                   , na.rm = T)}
means_vector <- as.vector(unlist(result))
analysis_df <- cbind(colony_buff_refined, Plastic_debris_mean = means_vector)

# Ordering colony-wise exposure in an intact order ----
analysis_df$Colony <- factor(analysis_df$Colony, levels = analysis_df$Colony)

# Plotting a basic bar plot for analysis_df ----
analysis_df_bar_plot <- ggplot(data = analysis_df, aes(x = Colony, y = Plastic_debris_mean)) +
  geom_col(colour = "black", fill = "skyblue") +
  scale_y_continuous(name = "Plastic debris mean") +
  coord_cartesian(ylim = c(1.1,5)) +
  theme_classic() +
  theme(axis.line = element_line(colour = "grey"), axis.text.x = element_text(angle = 90)) 
analysis_df_bar_plot

# Defining a standard deviation and standard error of mean column in analysis_df for error bars in plots ----
sd_result <- vector("list", 11)
for(i in 1:11){sd_result[[i]] <- sd(values
                                    (
                                      crop
                                      (plastics, colony_buff_refined[i,2])
                                    )
                                    , na.rm = T)}
sd_vector <- as.vector(unlist(sd_result))
analysis_df_2 <- cbind(analysis_df, Standard_deviation = sd_vector)\
View(analysis_df_2)

# Defining a sample size column in analysis_df_2 ----
ss_result <- vector("list", 11)
for(i in 1:11){
  ss_result[[i]] <- sum(!is.na
                        (values(crop(plastics, colony_buff_refined[i,2])))
  )
}
Sample_size <- as.vector(unlist(ss_result))
analysis_df_2 <- cbind(analysis_df, Sample_size)
transmute_seom <- analysis_df_2 |> transmute(Standard_error_of_mean = Standard_deviation/sqrt(Sample_size))
analysis_df_2 <- cbind(analysis_df_2, transmute_seom)


# Plot with standard deviation ----
analysis_df_box_plot <- ggplot(data = analysis_df_2, aes(x = Colony, y = Plastic_debris_mean,
                                                         ymin = Plastic_debris_mean - Standard_deviation,
                                                         ymax = Plastic_debris_mean + Standard_deviation)) +
  scale_y_continuous(name = "Plastic debris mean") +
  geom_boxplot() +
  theme_classic() +
  theme(axis.line = element_line(colour = "grey"), axis.text.x = element_text(angle = 90)) +
  geom_errorbar()
analysis_df_box_plot

# Plot with seom ----
analysis_df_seom <- analysis_df_box_plot <- ggplot(data = analysis_df_2, aes(x = Colony, y = Plastic_debris_mean,
                                                                             ymin = Plastic_debris_mean - Standard_error_of_mean,
                                                                             ymax = Plastic_debris_mean + Standard_error_of_mean)) +
  scale_y_continuous(name = "Plastic debris mean") +
  geom_boxplot() +
  theme_classic() +
  theme(axis.line = element_line(colour = "grey"), axis.text.x = element_text(angle = 90)) +
  geom_errorbar()
analysis_df_seom


# Adding a sample size column in analysis_df ----
ss_result <- vector("list", 11)
for(i in 1:11){
  ss_result[[i]] <- sum(!is.na
                        (values(crop(plastics, colony_buff_refined[i,2])))
  )
}
Sample_size <- as.vector(unlist(ss_result))
analysis_df_2 <- cbind(analysis_df, Sample_size)
View(analysis_df_2)

# Normality testing: ----

# 1. Density plots ----
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

# 2. QQ plots ----
qplot <- vector('list', 11)
for(i in 1:11){ 
  qplot[[i]] <- local({
    i <- i
    print(ggqqplot(na.omit(values(crop(plastics, colony_buff_refined[i,2])))), 
          main = "Density plot",
          )
    }
    )
}
library(gridExtra)
for(j in 1:11){
  grid.arrange(grobs = qplot, ncol = 5, nrow = 3)
}

# 3. Shapiro-Wilk's test ----
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



# Creating a boxplot for easy visualization ----

boxplot_plastic_debris_values <- ggplot(data = new_analysis_df, aes(x = Colony, y = Plastic_debris_value)) +
  geom_point() +
  geom_boxplot() +
  theme_classic() +
  scale_y_continuous(name = "Plastic debris value") +
  theme(axis.line = element_line(colour = "grey"), axis.text.x = element_text(angle = 90)) 
boxplot_plastic_debris_values

# Playing around with linear models----
# Checking the assumptions of lm with categorical predictors 

install.packages("performance")
library(performance)
install.packages("see")
library(see)
trial_lm <- lm(data = new_analysis_df, Plastic_debris_value ~ Colony)
trial_lm
summary(trial_lm)
anova(trial_lm)
plot(trial_lm)
check_model(trial_lm)
shapiro.test(trial_lm$residuals) # Residuals are non-normal
library(dplyr)
trimmed_new_analysis_df <- new_analysis_df |> filter(!grepl('Alkefjellet', Colony))
trimmed_lm <- lm(data = trimmed_new_analysis_df, Plastic_debris_value ~ Colony)
shapiro.test(trimmed_lm$residuals) # Residuals are non-normal

