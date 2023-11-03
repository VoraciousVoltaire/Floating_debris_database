library(stars)
library(ggplot2)

# LOADING RASTER
dtm_harv <- read_stars("data/HARV/HARV_dtmCrop.tif")

# PLOTTING
ggplot () +
  geom_stars(data = dtm_harv) +
  scale_fill_viridis_c() # c stands for continuous, d for discrete; viridis takes care of colour-blindness

# DSM is Digital Surface Model, DTM is Digital Terraain Model, CHM is Canopy Height Model
# DSM - DTM = CHM

# LOADING VECTOR
library(sf)
plots_harv <- st_read("data/HARV/harv_plots.shp")
ggplot() +
  geom_sf(data = plots_harv, mapping = aes(color = plot_type))

# Use st_crs for looking up coordinate reference systems of vector or raster objects
# UTM stands for Universal Transverse Mercator is most commonly used for ecological research 

plots_harv_utm <- st_transform(dtm_harv, st_crs(plots_harv))
plot_elevations <- aggregate(dtm_harv, plots_harv_utm, mean, as_points = F)

# Use facet_wrap function to form different sub-plots and use the tilda sign inside it to specify the feature on which basis the sub-plots should be designed upon.
# By default, geom_sf converts units to lats and longs. To counter this, use coord_sf in ggplot to set the projection system under the command datum = 

# To convert geolocator data into a shapefile,
harv_plots <- st_read("data/HARV/harv_plots.csv",
                      options= c(X_POSSIBLE_NAMES = longitude,
                                 Y_POSSIBLE_NAMES = latitude),
                      crs = 4326)
st_write(harv_plots, "harv_plots_new.shp")

# To crop raster data to constrain it only to the limits of the overlying vector data:
harv_boundary <- st_read("data/HARV/harv_boundary.shp")
dtm_harv_cropped <- st_crop(dtm_harv, harv_boundary) # if we pass crop = F as one of the arguments in st_crop, the original raster dimensions would be preserved but masking instead of cropping would take place
ggplot() +
  geom_stars(data = dtm_harv_cropped) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_sf(data = harv_boundary, alpha = 0)

bbox = st_bbox(c(xmin=731000,xmax=732000,ymin=4713000,ymax=4714000), st_crs(dtm_harv))
harv_dtm_small <- st_crop(dtm_harv, bbox)

write_stars(harv_dtm_small, "harv_dtm_small.tif")
read_stars("harv_dtm_small.tif")

# A shapefile is made up of a bunch of files with extensions such as .shp, .shx, .dbf and .prj